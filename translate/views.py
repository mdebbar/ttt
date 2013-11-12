from django.shortcuts import render


def index(req):
  context = {
    'test': 'TEST'
  }
  return render(req, 'translate/index.html', context)
