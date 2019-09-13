from django.shortcuts import render
from rest_framework import generics

from tool.models import BackgroundFormat, Job
from tool.api.serializers import BackgroundFormatSerializer, JobSerializer

class BackgroundFormatListCreate(generics.ListCreateAPIView):
	queryset = BackgroundFormat.objects.all()
	serializer_class = BackgroundFormatSerializer

class JobListCreate(generics.ListCreateAPIView):
	queryset = Job.objects.all()
	serializer_class = JobSerializer