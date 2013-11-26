from django.conf.urls import patterns, include, url

from django.shortcuts import render

def index(req):
  return render(req, 'index.html')

urlpatterns = patterns('',
    url(r'^$', index, name='index'),
    url(r'^blast/', include('blast.urls')),
)
