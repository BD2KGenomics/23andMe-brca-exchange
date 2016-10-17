# 23andMe-brca-exchange
Get expert, curated information about your 23andMe results in the BRCA1 and 2 regions.

# Instructions

When you sign up to be a developer with 23andMe, you'll get a `client_id`, and a `client secret` (keep the secret safe).

You can start the server by supplying your client ID:

```
python client.py --client-id=5f6d821a3d48c403fd931417f5711ad9
```

You'll then be prompted for your `client_secret`, which you should enter.

If all goes well, the server will be listening on http://localhost:5000, and you can go ahead and visit it.

This repository has some [modifications](https://github.com/23andMe/api-example-flask) to include the official 23andMe [branded button](https://api.23andme.com/docs/jslib).
