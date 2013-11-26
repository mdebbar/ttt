from django.conf.urls import patterns, include, url

from django.shortcuts import render
from blast import api


def index(req):
  algos = sorted(api.DATABASES.iterkeys())
  return render(req, 'index.html', {'algos': algos})

urlpatterns = patterns('',
    url(r'^$', index, name='index'),
    url(r'^blast/', include('blast.urls')),
)
