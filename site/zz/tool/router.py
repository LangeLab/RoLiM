from tool.api import viewsets
from rest_framework.routers import DefaultRouter

router = DefaultRouter()
router.register(r'job', JobViewSet, basename='job')
urlpatterns = router.urls