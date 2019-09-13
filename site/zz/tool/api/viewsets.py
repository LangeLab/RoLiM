import uuid

#from django.views.decorators.csrf import ensure_csrf_cookie
from rest_framework import viewsets
from rest_framework.exceptions import ParseError
from rest_framework.parsers import FileUploadParser, MultiPartParser, FormParser
from rest_framework.response import Response
from rq import Queue
from redis import Redis
import django_rq

from tool.models import (Job, DataFormat, BackgroundFormat)
from tool.api.serializers import (JobSerializer, DataFormatSerializer,
									BackgroundFormatSerializer)
from tool.tasks import new_job



#@ensure_csrf_cookie
class JobViewSet(viewsets.ModelViewSet):
	serializer_class = JobSerializer
	parser_class = (MultiPartParser, FormParser,)

	def get_queryset(self):
		# Restrict queryset to jobs associated with user session.
		queryset = Job.objects.filter(session_id=self.request.session.session_key)
		return queryset

	def perform_create(self, serializer):
		while True:
			self.jobcode = uuid.uuid4().hex
			if not Job.objects.filter(jobcode=self.jobcode).exists():
				break
		
		try:
			self.request.session['jobcodes'].append(self.jobcode)
		except KeyError:
			self.request.session['jobcodes'] = [self.jobcode]

		serializer.save(session_id=self.request.session.session_key,
						jobcode=self.jobcode)
		django_rq.enqueue(new_job, self.jobcode)


class BackgroundFormatViewSet(viewsets.ModelViewSet):
	queryset = BackgroundFormat.objects.all()
	serializer_class = BackgroundFormatSerializer