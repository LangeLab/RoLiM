import uuid

from rest_framework import viewsets
from rest_framework.exceptions import ParseError
from rest_framework.parsers import FileUploadParser, MultiPartParser, FormParser
from rq import Queue
from redis import Redis
import django_rq

from patterndetection.tool.models import (
	Job, ForegroundFormat, ContextFormat, ExtensionDirection,
	RedundancyLevel, OriginalRowMerge
)
from patterndetection.tool.serializers import (
	JobSerializer,
	ForegroundFormatSerializer,
	ContextFormatSerializer,
	ExtensionDirectionSerializer,
	RedundancyLevelSerializer,
	OriginalRowMergeSerializer,
)
from patterndetection.tool.tasks import new_job


class JobViewSet(viewsets.ModelViewSet):
	"""Views for Job model REST endpoints."""
	serializer_class = JobSerializer
	parser_class = (MultiPartParser, FormParser)

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

		# Save new job to database.
		serializer.save(
			session_id=self.request.session.session_key,
			jobcode=self.jobcode
		)

		# Add job to task queue.
		queue = django_rq.get_queue('default')
		queue.enqueue(new_job, args=(self.jobcode,), job_timeout='48h')


class ForegroundFormatViewSet(viewsets.ModelViewSet):
	"""Views for ForegroundFormat model REST endpoints."""
	queryset = ForegroundFormat.objects.all()
	serializer_class = ForegroundFormatSerializer


class ContextFormatViewSet(viewsets.ModelViewSet):
	"""Views for ForegroundFormat model REST endpoints."""
	queryset = ContextFormat.objects.all()
	serializer_class = ContextFormatSerializer


class ExtensionDirectionViewSet(viewsets.ModelViewSet):
	"""Views for ExtensionDirection model REST endpoints."""
	queryset = ExtensionDirection.objects.all()
	serializer_class = ExtensionDirectionSerializer


class RedundancyLevelViewSet(viewsets.ModelViewSet):
	"""Views for ExtensionDirection model REST endpoints."""
	queryset = RedundancyLevel.objects.all()
	serializer_class = RedundancyLevelSerializer


class OriginalRowMergeViewSet(viewsets.ModelViewSet):
	"""Views for ExtensionDirection model REST endpoints."""
	queryset = OriginalRowMerge.objects.all()
	serializer_class = OriginalRowMergeSerializer