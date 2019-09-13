from django.contrib import admin
from django.urls import path, include
from django.conf.urls.static import static
from django.conf import settings
#from .router import router

urlpatterns = [
    path('admin/', admin.site.urls),
    #path('api/', include(router.urls)),
    path('', include('frontend.urls')),
    path('', include('tool.urls')),
]

urlpatterns += [
    path('django-rq/', include('django_rq.urls'))
]

if settings.DEBUG:
	urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)