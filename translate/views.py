from django.shortcuts import render


def index(req):
  return render(req, 'translate/index.html')
