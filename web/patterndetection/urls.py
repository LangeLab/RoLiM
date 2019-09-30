from django.contrib import admin
from django.urls import path, include
from django.conf.urls.static import static
from django.conf import settings
from rest_framework import routers
from patterndetection.tool import viewsets

router = routers.DefaultRouter()
router.register(r'job', viewsets.JobViewSet, basename='job')
router.register(r'foregroundformat', viewsets.ForegroundFormatViewSet, basename='foregroundformat')
router.register(r'contextformat', viewsets.ContextFormatViewSet, basename='contextformat')
router.register(r'extensiondirection', viewsets.ExtensionDirectionViewSet, basename='extensiondirection')

urlpatterns = [
    path('admin/', admin.site.urls),
    path('api/', include(router.urls)),
    path('', include('patterndetection.frontend.urls')),
    path('api-auth/', include('rest_framework.urls', namespace='rest_framework'))
]
