from django.urls import path
from django.views.generic import RedirectView
from . import views

urlpatterns = [
	path('rolim/', views.index),
	path('rolim/textfile/', views.prealigned_text_file),
	path('rolim/peptidelist/', views.text_file_peptide_list),
	path('rolim/compoundresiduefile/', views.compound_residue_file),
	path('', RedirectView.as_view(url='https://langelab.med.ubc.ca/')),
]