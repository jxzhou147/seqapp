from django.db import models
from .functions import sequence_len


class Sequence(models.Model):
    sequence_text = models.CharField(max_length=200)
    pub_date = models.DateTimeField('date published')

    def result(self):
        return sequence_len(self.sequence_text)

    def __str__(self):
        return self.sequence_text
