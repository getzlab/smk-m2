#!/bin/bash

auth="$HOME"/"auth.json"

if [ -f $auth ]; then
    echo "Authorized"
else
# auth
    curl --request POST \
        --url https://us-central1-broad-dsde-prod.cloudfunctions.net/martha_v2 \
        --header "authorization: Bearer $(gcloud auth print-access-token)" \
        --header 'content-type: application/json' \
        --data "{ \"url\": \"$1\" }" | /demo-mount/tools/jq-linux64 '.googleServiceAccount.data' > $auth
fi

