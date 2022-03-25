from django.db import models


class Function(models.Model):
    function_name = models.CharField(max_length=200)
    function_path = models.FileField()
    pub_date = models.DateTimeField('data published')

    def __str__(self):
        return self.function_name
