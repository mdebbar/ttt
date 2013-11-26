from django.conf.urls import patterns, url

from blast import views

urlpatterns = patterns('',
  url(r'^$', views.blast, name='blast'),
)
