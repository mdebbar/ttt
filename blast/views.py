import json, api
from django.http.response import HttpResponse


# responds to ajax requests
def blast(req):
  #response = api.blast(req.POST['seq'], req.POST.get('algo'))
  response = api.sample()
  return HttpResponse(json.dumps(response))
