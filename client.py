import getpass
import logging
import os
from optparse import OptionParser

import requests
import flask
from flask import request
from requests_oauthlib import OAuth2Session

import ga4gh.client as g4client


PORT = 5000
API_SERVER = 'api.23andme.com'
BASE_CLIENT_URL = 'http://localhost:%s/' % PORT
DEFAULT_REDIRECT_URI = '%sapp/' % BASE_CLIENT_URL
PAGE_HEADER = "23andMe + GA4GH"

# So we don't get errors if the redirect uri is not https.
os.environ["OAUTHLIB_INSECURE_TRANSPORT"] = '1'

# Pass in more scopes through the command line, or change these.
DEFAULT_SNPS = ['rs12913832', 'rs3088053', 'rs1000068']
DEFAULT_SCOPES = ['names', 'basic'] + DEFAULT_SNPS

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
parser.add_option("-r", "--redirect_uri", dest="redirect_uri", default=DEFAULT_REDIRECT_URI,
                  help="Your client's redirect_uri [%s]" % DEFAULT_REDIRECT_URI)
parser.add_option("-a", "--23andMe-api-server", dest="t23andMe_api_server", default=API_SERVER,
                  help="Almost always: [api.23andme.com]")
parser.add_option("-p", "--select-profile", dest='select_profile', action='store_true', default=False,
                  help='If present, the auth screen will show a profile select screen')
parser.add_option("-f", "--ga4gh-api-server", dest="ga4gh_api_server", help="The GA4GH API server location.")
parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False,
                  help="Whether or not to provide debugging output.")

(options, args) = parser.parse_args()

DEBUG = options.debug

BASE_API_URL = "https://%s" % options.t23andMe_api_server
API_AUTH_URL = '%s/authorize' % BASE_API_URL
API_TOKEN_URL = '%s/token/' % BASE_API_URL

API_SERVER_GA4GH = options.ga4gh_api_server

if options.select_profile:
    API_AUTH_URL += '?select_profile=true'

redirect_uri = options.redirect_uri
client_id = options.client_id
if options.client_secret:
    client_secret = options.client_secret

scopes = options.scopes or DEFAULT_SCOPES

if not options.client_id:
    print "missing param: CLIENT_ID:"
    parser.print_usage()
    print "Please navigate to your developer dashboard [%s/dev/] to retrieve your client_id.\n" % BASE_API_URL
    exit()

if not client_secret:
    print "Please navigate to your developer dashboard [%s/dev/] to retrieve your client_secret." % BASE_API_URL
    client_secret = getpass.getpass("Please enter your client_secret: ")

app = flask.Flask(__name__)


@app.route('/')
def index():
    """Here, we authenticate the user before transitioning to the app.  There
    should be no way of getting to the app without this step."""
    ttam_oauth = OAuth2Session(client_id, redirect_uri=redirect_uri,
                               scope=scopes)
    auth_url, state = ttam_oauth.authorization_url(API_AUTH_URL)
    print("Authentication URL: %s" % auth_url)
    print("State: %s" % state)
    demo_url = "http://localhost:5000/demo/"
    return flask.render_template('index.html', auth_url=auth_url,
        page_header=PAGE_HEADER, page_title=PAGE_HEADER, client_id=client_id,
        demo_url=demo_url)


@app.route('/demo/')
def demo():
    """This is the view responsible for handling the case where the user wants
    to demo the application, without logging in to 23andMe."""
    genotype_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/genotype/"),
                                     params={'locations': ' '.join(DEFAULT_SNPS)},
                                     verify=False)
    basic_response = requests.get("%s%s" % (BASE_API_URL, "/1/demo/user/"),
                                     verify=False)
    if basic_response.status_code == 200:
        print(basic_response)
    else:
        basic_response.raise_for_status()

    if genotype_response.status_code == 200:
        return flask.render_template('app.html', page_header=PAGE_HEADER,
            response_json=genotype_response.json(), home_url=BASE_CLIENT_URL,
            page_title=PAGE_HEADER, client_id=client_id)
    else:
        genotype_response.raise_for_status()


def _23andMe_queries(client_id, client_secret, redirect_uri):
    """Handles interaction with the 23andMe API.  Returns the data."""
    # Hit the /token endpoint to get the access_token.
    ttam_oauth = OAuth2Session(client_id, redirect_uri=redirect_uri)
    token_dict = ttam_oauth.fetch_token(API_TOKEN_URL,
                                        client_secret=client_secret,
                                        authorization_response=request.url)

    access_token = token_dict['access_token']

    headers = {'Authorization': 'Bearer %s' % access_token}

    genotype_response = requests.get("%s%s" % (BASE_API_URL, "/1/genotype/"),
                                    params={'locations': ' '.join(DEFAULT_SNPS)},
                                    headers=headers,
                                    verify=False)
    basic_response = requests.get("%s%s" % (BASE_API_URL, "/1/user/"),
                                    headers=headers,
                                    verify=False)
    #names_response = requests.get("%s%s" % (BASE_API_URL, "/1/names/"),
    #                                 headers=headers,
    #                                 verify=False)
    #if names_response.status_code == 200:
    #    print(names_response)
    #else:
    #    names_response.raise_for_status()
    return genotype_response, basic_response

def _ga4gh_queries():
    """Performs queries against the GA4GH server."""
    if DEBUG:
        httpClient = g4client.HttpClient(API_SERVER_GA4GH, logLevel=logging.DEBUG)
    else:
        httpClient = g4client.HttpClient(API_SERVER_GA4GH)
    datasets = list(httpClient.search_datasets())
    variant_sets = list(httpClient.search_variant_sets(dataset_id=datasets[0].id))
    iterator = httpClient.search_variants(variant_set_id=variant_sets[0].id,
        reference_name="1", start=45000, end=50000)
    results = set()
    for variant in iterator:
        r = (variant.reference_name, variant.start, variant.end,\
            variant.reference_bases, variant.alternate_bases)
        results.add(r)
    return results

@app.route('/app/')
def app2():
    """Represents our application, which makes use of 2 APIs: 23andMe, and
    BRCA Exchange (via GA4GH)."""
    # Query the 2 APIs and get data responses.
    genotype_response, basic_response = _23andMe_queries(client_id, client_secret, redirect_uri)
    results = _ga4gh_queries()

    for v in results:
        print(v)

    if basic_response.status_code == 200:
        print(basic_response)
    else:
        basic_response.raise_for_status()
    if genotype_response.status_code == 200:
        if 'code' in request.args.to_dict():
            code = request.args.to_dict()['code']
        else:
            code = None
        return flask.render_template('app.html', page_header=PAGE_HEADER,
            response_json=genotype_response.json(), home_url=BASE_CLIENT_URL,
            page_title=PAGE_HEADER, client_id=client_id, code=code,
            ga4gh_results=results)
    else:
        genotype_response.raise_for_status()


if __name__ == '__main__':
    print "A local client for the Personal Genome API is now initialized."
    app.run(debug=DEBUG, port=PORT)
