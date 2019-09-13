from django.shortcuts import render
from rest_framework import viewsets
from rest_framework.parsers import MultiPartParser, FormParser
from .models import DataUpload, BackgroundUpload
from .serializers import DataUploadSerializer

class DataUploadViewSet(viewsets.ModelViewSet):
	"""ViewSet for uploaded files."""

	queryset = DataUpload.objects.all()
	serializer_class = DataUploadSerializer
	parser_classes = (MultiPartParser, FormParser,)

	def perform_create(self, serializer):
		serializer.save(file_path=self.request.data.get('file_path'))

class BackgroundUploadViewSet(viewsets.ModelViewSet):
	pass

class DataTypesViewSet(viewsets.ModelViewSet):
	pass

class BackgroundTypesViewSet(viewsets.ModelViewSet):
	pass

class ParametersViewSet(viewsets.ModelViewSet):
	pass

class MetaDataViewSet(viewsets.ModelViewSet):
	pass

class ProteomeViewSet(viewsets.ModelViewSet):
	pass

class JobViewSet(viewsets.ModelViewSet):
	pass
