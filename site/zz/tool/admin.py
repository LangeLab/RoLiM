from django.contrib import admin
from tool.models import (Job, Pattern, Figure, DataFormat, BackgroundFormat,)

# Register your models here.
admin.site.register(Job)
admin.site.register(Pattern)
admin.site.register(Figure)
admin.site.register(DataFormat)
admin.site.register(BackgroundFormat)