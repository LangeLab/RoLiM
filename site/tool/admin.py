from django.contrib import admin
from tool.models import (DataUpload, BackgroundUpload, MetaData,
							Proteome, DataTypes, Parameters, Job)


# Register your models here.
admin.site.register(MetaData)
admin.site.register(Proteome)
admin.site.register(DataTypes)
admin.site.register(Parameters)
admin.site.register(DataUpload)
admin.site.register(BackgroundUpload)
admin.site.register(Job)