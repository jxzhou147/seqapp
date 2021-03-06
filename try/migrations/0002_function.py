# Generated by Django 3.2.5 on 2022-03-24 15:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('try', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Function',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('function_name', models.CharField(max_length=200)),
                ('function_path', models.FileField(upload_to='')),
                ('pub_date', models.DateTimeField(verbose_name='data published')),
            ],
        ),
    ]
