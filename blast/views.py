import json, api
from django.http.response import HttpResponse


# responds to ajax requests
def blast(req):
  #response = api.sample()
  response = api.blast(req.POST['seq'], req.POST.get('algo'))
  return HttpResponse(json.dumps(response))
