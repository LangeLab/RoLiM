from django.urls import path
from django.conf.urls.static import static
from django.conf import settings
from django.views.generic.base import RedirectView

from . import views

urlpatterns = [
	path('', views.index),
	#path('', RedirectView.as_view(url='static/index.html', permanent=False)),
]