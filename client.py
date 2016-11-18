"""
What was done,
what we plan to do
TODOs
"""


import getpass
import logging
import sys
import re
import pickle
import json
import csv
import os
from optparse import OptionParser

import google.protobuf.json_format as json_format
import requests
import flask
from flask import request
from requests_oauthlib import OAuth2Session

from ga4gh_client import client as g4client
from ga4gh.exceptions import RequestNonSuccessException


PORT = 5000
API_SERVER = 'api.23andme.com'
BASE_CLIENT_URL = 'http://localhost:%s/' % PORT
DEFAULT_APP_REDIRECT_URI = '%sapp/' % BASE_CLIENT_URL
DEFAULT_API_REDIRECT_URI = '%svariants/search/' % BASE_CLIENT_URL
PAGE_HEADER = "23andMe + GA4GH"
#REFERENCE_NAMES = [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']
REFERENCE_NAMES = ["13", "17"]
BRCA2_START = 32889611
access_token = None

# So we don't get errors if the redirect uri is not https.
os.environ["OAUTHLIB_INSECURE_TRANSPORT"] = '1'

# Pass in more scopes through the command line, or change these.
DEFAULT_SNPS = ['rs12913832', 'rs3088053', 'rs1000068', 'rs206118', 'rs206115', 'rs3094315', 'i3000001']
DEFAULT_SCOPES = ['names', 'basic', 'email', 'analyses', 'genomes'] + DEFAULT_SNPS

# The program will ask for a client_secret if you choose to not hardcode one
# here.
client_secret = None

parser = OptionParser(usage="usage: %prog -i CLIENT_ID [options]")
parser.add_option("-i", "--client-id", dest="client_id", default='',
                  help="Your client_id [REQUIRED]")
parser.add_option('-s', '--scopes', dest='scopes', action='append', default=[],
                  help='Your requested scopes. Eg: -s basic -s rs12913832')
parser.add_option("-c", "--client-secret", dest='client_secret',
                  help='The client secret')
parser.add_option("-r", "--app-redirect-uri", dest="app_redirect_uri", default=DEFAULT_APP_REDIRECT_URI, help="Your client's redirect_uri [%s]" % DEFAULT_APP_REDIRECT_URI)
parser.add_option("-z", "--api-redirect-uri", dest="api_redirect_uri", default=DEFAULT_API_REDIRECT_URI, help="Your client's redirect_uri [%s]" % DEFAULT_API_REDIRECT_URI)
parser.add_option("-a", "--23andMe-api-server", dest="t23andMe_api_server", default=API_SERVER, help="Almost always: [api.23andme.com]")
parser.add_option("-p", "--select-profile", dest='select_profile', action='store_true', default=False, help='If present, the auth screen will show a profile select screen')
parser.add_option("-f", "--ga4gh-api-server", dest="ga4gh_api_server", help="The GA4GH API server location.")
parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False,
                  help="Whether or not to provide debugging output.")
parser.add_option("-k", "--snps-data", dest='snps_data',
                  help='A SNPS data file to use.')

(options, args) = parser.parse_args()


def _23andMe_queries(client_id, client_secret, app_redirect_uri):
    """Handles interaction with the 23andMe API.  Returns the data."""
    global access_token
    if not access_token:
        ttam_oauth = OAuth2Session(client_id, redirect_uri=app_redirect_uri)
        token_dict = ttam_oauth.fetch_token(API_TOKEN_URL,
                                            client_secret=client_secret,
                                            authorization_response=request.url)

        access_token = token_dict['access_token']

    headers = {'Authorization': 'Bearer %s' % access_token}

    # The documentation says you shouldn't call 'names' if the user wants to be
    # anonymous.  We should check on this.  The only reason why we're using it
    # now is to get the demo profiles.
    names_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/names/"),
                                    headers=headers,
                                    verify=True)

    genotype_responses = []
    for profile in names_response.json()['profiles']:
        genotype_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/genotypes/%s/" % profile['id']),
                                        params={'locations': ' '.join(locations), 'format': 'embedded'},
                                        headers=headers,
                                        verify=True)
        genotype_responses.append(genotype_response)
    user_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/user/"),
                                    headers=headers,
                                    verify=True)
    #profilepic_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/profile_picture/SP1_FATHER_V4/"),
    #                                headers=headers,
    #                                verify=True)
    #family_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/family_members/"),
    #                                headers=headers,
    #                                verify=True)
    return genotype_responses, user_response, names_response, None, None, None, None


def _g4_queries():
    """Performs queries against the GA4GH server."""
    if DEBUG:
        httpClient = g4client.HttpClient(API_SERVER_GA4GH, logLevel=logging.DEBUG)
    else:
        httpClient = g4client.HttpClient(API_SERVER_GA4GH)
    # There is currently only 1 dataset available in BRCA, but we'll be robust
    # and iterate as if there were more.
    datasets = list(httpClient.search_datasets())
    results = list()
    for dataset in datasets:
        # There should be 3 variant sets; we're only concerned with hg37 for
        # now though.
        variant_sets = list(httpClient.search_variant_sets(dataset_id=dataset.id))
        variant_set = filter(lambda x: x.id == 'brca-hg37', variant_sets)[0]
        for reference_name in REFERENCE_NAMES:
            iterator = httpClient.search_variants(variant_set_id=variant_set.id, reference_name=reference_name, start=BRCA2_START, end=BRCA2_START+2000)
            #iterator = httpClient.search_variants(variant_set_id=variant_set.id, reference_name=reference_name, start=105598600, end=105598700)
                #reference_name=reference_name, start=32315650, end=32315660)
                #reference_name="13", start=0, end=500000)
            for variant in iterator:
                r = (variant.reference_name, variant.start, variant.end,\
                    variant.reference_bases, variant.alternate_bases, variant.id,\
                    variant.info, variant.names)
                #r = variant
                results.append(r)
    return (datasets, variant_sets, results)


def _compute_locations_from_snps_file(start=41196311, end=41196314, reference_name="13", s=""):
    """Computes a more reasonable list of SNPs than DEFAULT_SNPS.

    It returns the rsIDs from the given SNPs file that are associated with the
    given reference name, and has a position that falls within the given
    range."""
    cross = []
    cross_augmented = []
    found_header = False
    with open(s, 'r') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if not row[0].startswith('#'):
                if not found_header:
                    found_header = True
                    continue
                index, snp, ch, p = row[0], row[1], row[2], int(row[3])
                if ch == reference_name and p > start and p < end:
                    cross.append(snp)
                    cross_augmented.append((snp, p))
    print "Crosses: %s" % len(cross)
    #return cross if len(cross) > 0 else ' '.join(DEFAULT_SNPS)
    #return cross + ' '.join(DEFAULT_SNPS) + scopes
    #return " ".join([x[0] for x in cross])
    return cross, cross_augmented


if not options.snps_data:
    print("Should specify --snps-data option.")
    sys.exit(1)
SNPS_DATA_FILE = options.snps_data

locations, locations_augmented = _compute_locations_from_snps_file(start=32889611, end=32889611+2000, s=SNPS_DATA_FILE)
#locations = _compute_locations_from_snps_file(start=41196311, end=41277500, s=SNPS_DATA_FILE)
#locations = _compute_locations_from_snps_file(start=41196311, end=41196314, s=SNPS_DATA_FILE)
print "Locations: %s %s" % (len(locations), locations)
scopes = options.scopes or (DEFAULT_SCOPES + locations)
print "Scopes: %s %s" % (len(scopes), scopes)

DEBUG = options.debug

BASE_API_URL = "https://%s" % options.t23andMe_api_server
API_AUTH_URL = '%s/authorize' % BASE_API_URL
API_TOKEN_URL = '%s/token/' % BASE_API_URL

API_SERVER_GA4GH = options.ga4gh_api_server

if options.select_profile:
    API_AUTH_URL += '?select_profile=true'

app_redirect_uri = options.app_redirect_uri
api_redirect_uri = options.api_redirect_uri
client_id = options.client_id
if options.client_secret:
    client_secret = options.client_secret

if not options.client_id:
    print "missing param: CLIENT_ID:"
    parser.print_usage()
    print "Please navigate to your developer dashboard [%s/dev/] to retrieve your client_id.\n" % BASE_API_URL
    exit()

if not client_secret:
    print "Please navigate to your developer dashboard [%s/dev/] to retrieve your client_secret." % BASE_API_URL
    client_secret = getpass.getpass("Please enter your client_secret: ")

app = flask.Flask(__name__)
app.config['SECRET_KEY'] = 'abc'    # May not need this here?  It wasn't there


@app.route('/')
def index():
    """Here, we authenticate the user before transitioning to the app.  There
    should be no way of getting to the app without this step."""
    #start = int(flask.request.args.get('start', 32889611))
    #end = int(flask.request.args.get('end', 32889611 + 1000))
    #reference_name = str(flask.request.args.get('reference_name', "13"))
    #ttam_oauth = OAuth2Session(client_id, redirect_uri="http://localhost:5000/app/?start={}&end={}&reference_name={}".format(start, end, reference_name),
    ttam_oauth = OAuth2Session(client_id, redirect_uri=app_redirect_uri, scope=scopes)
    auth_url, state = ttam_oauth.authorization_url(API_AUTH_URL)
    return flask.render_template('index.html', auth_url=auth_url,
        page_header=PAGE_HEADER, page_title=PAGE_HEADER, client_id=client_id)


@app.route('/variants/search/')
def variants_search_endpoint():
    flask.session['23code'] = flask.request.args.get('code')
    start = int(flask.request.args.get('start', 32889611))
    end = int(flask.request.args.get('end', 32889611 + 1000))
    reference_name = str(flask.request.args.get('reference_name', "13"))
    
    # Get GA4GH variants
    g4 = g4client.HttpClient(API_SERVER_GA4GH)
    # could take a long time for large regions beware, maybe kick out if a too
    # large region is requested
    variants = list(g4.search_variants(variant_set_id="brca-hg37", start=start,
        end=end, reference_name=reference_name))

    locations, _ = _compute_locations_from_snps_file(start=start, end=end,
        reference_name=reference_name, s=SNPS_DATA_FILE)

    print(locations)
    global access_token
    if not access_token:
        ttam_oauth = OAuth2Session(client_id, redirect_uri=api_redirect_uri)
        token_dict = ttam_oauth.fetch_token(API_TOKEN_URL,
            client_secret=client_secret, authorization_response=request.url)

        access_token = token_dict['access_token']
    
    headers = {'Authorization': 'Bearer %s' % access_token}

    genotype_responses = []
    for profile in names_response:
        genotype_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/genotypes/%s/" % profile['id']),
                                        params={'locations': locations, 'format': 'embedded'},
                                        headers=headers,
                                        verify=True)
        genotype_responses.append(genotype_response)

    return flask.jsonify({"23andme": genotype_responses.json(), "g4": [json_format._MessageToJsonObject(v, True) for v in variants], "locations": locations, "query": {"start": start, "end": end, "reference_name": reference_name}})


@app.route('/app/')
def app2():
    """Represents our application, which makes use of 2 APIs: 23andMe, and BRCA
    Exchange (via GA4GH)."""
    genotype_responses, user_response, names_response, profilepic_response, family_response, neanderthal_response, relatives_response = _23andMe_queries(client_id, client_secret, app_redirect_uri)

    # Compute some context data for later rendering.
    user_request_success = False
    if user_response.status_code == 200:
        user_request_success = True
    names_request_success = False
    if names_response.status_code == 200:
        names_request_success = True
    #profilepic_request_success = False
    #if profilepic_response.status_code == 200:
    #    profilepic_request_success = True
    #family_request_success = False
    #if family_response.status_code == 200:
    #    family_request_success = True
    #neanderthal_request_success = False
    #if neanderthal_response.status_code == 200:
    #    neanderthal_request_success = True
    #relatives_request_success = False
    #if relatives_response.status_code == 200:
    #    relatives_request_success = True
    if 'first_name' in names_response.json():
        account_first_name = names_response.json()['first_name']
    else:
        account_first_name = "first"
    if 'last_name' in names_response.json():
        account_last_name = names_response.json()['last_name']
    else:
        account_last_name = "last"
    if 'code' in request.args.to_dict():
        code = request.args.to_dict()['code']
    else:
        code = None

    # The algorithm that computes useful results with G4 and 23andMe data.
    temp = []
    g4 = g4client.HttpClient(API_SERVER_GA4GH)
    for reference_name in REFERENCE_NAMES:
        variants = list(g4.search_variants(variant_set_id="brca-hg37", start=BRCA2_START,
            end=BRCA2_START+1000, reference_name=reference_name))
        for gr in genotype_responses:
            profile = gr.json()
            r = None
            if 'genotypes' in profile:
                for call in profile['genotypes']:
                    for location in locations_augmented:
                        if call['location'] == location[0]:
                            for variant in variants:
                                if variant.start == location[1]:
                                    r = ("brca and 23andme have {}".format(location[0]),
                                        variant.info["Allele_Frequency"],
                                        "Individual presented: " + call['call'])
            if gr.status_code == 200:
                temp.append((True, gr.json(), r))
            else:
                temp.append((False, gr.json(), r))
    genotype_responses = temp

    # Reformat the genotype_responses list to include HTTP success status code;
    # and also, turn the second element into JSON (instead of a 'Request'
    # object).

    return flask.render_template('app.html', page_header=PAGE_HEADER,
        genotype_responses=genotype_responses,
        home_url=BASE_CLIENT_URL, user_response_json=user_response.json(),
        names_response_json=names_response.json(), page_title=PAGE_HEADER,
        client_id=client_id, code=code,
        user_request_success=user_request_success,
        names_request_success=names_request_success,
        api_results_url="http://localhost:5000/variants/search",
        #family_request_success=family_request_success,
        #neanderthal_request_success=neanderthal_request_success,
        #relatives_request_success=relatives_request_success,
        #profilepic_request_success=profilepic_request_success,
        account_first_name=account_first_name,
        account_last_name=account_last_name)


def _format_g4results(g):
    o = []
    for r in g:
        r = (r.id, r.names, r.reference_bases, r.reference_name, r.start, r.end, r.calls, r.info['Hg37_Start'].values[0].number_value, r.info['Hg37_End'].values[0].number_value, r.info['AFR_Allele_frequency_1000_Genomes'].values[0].string_value, r.info['EUR_Allele_frequency_1000_Genomes'].values[0].string_value, r.info['Chr'].values[0].string_value, r.info['Pathogenicity_expert'].values[0].string_value, r.info['Ref'].values[0].string_value, r.info['Alt'].values[0].string_value, r.info['Pos'].values[0].string_value, r.info['Allele_Frequency'].values[0].string_value, r.info['Gene_Symbol'].values[0].string_value)
        o.append(r)
    return o
#g4results = _format_g4results(g4results)


if __name__ == '__main__':
    app.run(debug=DEBUG, port=PORT)
