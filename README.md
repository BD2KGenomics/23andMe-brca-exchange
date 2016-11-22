# 23andMe-brca-exchange

An example web application for combining the results of 23andMe with that of BRCA Exchange.

23andMe is used to provide personal genomic results, and this service cross-references the information with that of BRCA Exchange.  You can use this to, say, tell if you're at particular risk of breast cancer.

# Using the service

The service wants/needs your permission to obtain information from 23andMe, in order to perform computations with it along with data from BRCA Exchange.  So before proceeding any further, you must log in to the 23andMe web site by following the given link, and you'll be redirected back to the application so you can see the results.

# Developer instructions

When you sign up to be a developer with 23andMe, you'll get a `client_id`, and a `client_secret` (keep the secret safe).

You can start the server by supplying your client ID:

```
python client.py --client-id 5f6d821a3d48c403fd931417f5711ad9
```

You'll then be prompted for your `client_secret`, which you should enter.

If all goes well, the server will be listening on http://localhost:5000, and you can go ahead and visit it with a web browser.


# About

This work was done by folks at the Genomics Institute, of University of California at Santa Cruz (UCSC).
