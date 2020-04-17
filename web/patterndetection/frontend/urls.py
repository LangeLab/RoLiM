from django.urls import path
from django.views.generic import RedirectView
from . import views

urlpatterns = [
	path('patterndetection/', views.index),
	path('patterndetection/textfile/', views.prealigned_text_file),
	path('patterndetection/peptidelist/', views.text_file_peptide_list),
	path('', RedirectView.as_view(url='https://langelab.med.ubc.ca/')),
]
