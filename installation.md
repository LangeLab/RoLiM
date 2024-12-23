## Python

RoLiM requires Python version 3.7.4 or later, however the latest stable release should be used unless prohibited by specific circumstances.

Instructions for installing Python from source can be found here: https://phoenixnap.com/kb/how-to-install-python-3-ubuntu

### Virtualenv and virtualenvwrapper

Install the virtualenv and virtualenvwrapper libraries for Python by following these instructions: https://itnext.io/virtualenv-with-virtualenvwrapper-on-ubuntu-18-04-goran-aviani-d7b712d906d5

Note that these resources should be installed system-wide, which will allow any version of Python installed on the system to be selected when creating a new virtual environment.

After installing, create a new virtual environment called "patterndetection" by running "mkvirtualenv patterndetection --python={path to preferred python executable}"

### RoLiM dependencies

Navigate to "~/project/RoLiM/" with "patterndetection" activated and install RoLiM dependencies using pip by running "python -m pip install -r requirements.txt"

### Installation notes

Python should be installed as an alt install rather than replacing the default system Python version in order to prevent possible system compatibility issues.

## MariaDB

RoLiM uses the MariaDB implementation of MySQL. Installation instructions can be found here: https://www.digitalocean.com/community/tutorials/how-to-install-mariadb-on-ubuntu-20-04

### Configuration

After installing, start a MySQL session by running ```"sudo mysql"```.

Next, create a database called "patterndetection" by running the following command:

```CREATE DATABASE IF NOT EXISTS patterndetection;```

Create a new user named "patterndetection" by running the following command:

```CREATE USER "patterndetection"@"localhost" identified by "{insert password from db_password file here}";```

Now grant this user full privileges on the new patterndetection database by running the following command:

```GRANT ALL PRIVILEGES ON patterndetection.* TO 'patterndetection'@'localhost';```

Create a second database called "merops" by running the following command:

```CREATE DATABASE IF NOT EXISTS merops;```

Now create a new user "zigzag" with full privileges on the merops database by running the following commands:

```
CREATE USER "zigzag"@"localhost" identified by "zigzagpass";
GRANT ALL PRIVILEGES ON merops.* TO 'patterndetection'@'localhost';
```

Finally, make sure that the new user privileges take effect by running the following command:

```FLUSH PRIVILEGES;```

### Installation notes

Python's MySQL connector requires the C client library API which can be installed by running "sudo apt-get install libmysqlclient-dev"

## Redis

Install Redis by following the instructions found at this address: https://www.digitalocean.com/community/tutorials/how-to-install-and-secure-redis-on-ubuntu-20-04

Open the Redis configuration file located at '/etc/redis/redis.conf' using a text editor. Find the line that begins with "require pass." Uncomment this line and replace the placeholder password in the file with the password in the redis_password file.

Confirm that a file exists at /etc/systemd/system/redis.service containing the following lines:

```
[Unit]
Description=Advanced key-value store
After=network.target
Documentation=http://redis.io/documentation, man:redis-server(1)

[Service]
Type=forking
ExecStart=/usr/bin/redis-server /etc/redis/redis.conf
PIDFile=/run/redis/redis-server.pid
TimeoutStopSec=0
Restart=always
User=redis
Group=redis
RuntimeDirectory=redis
RuntimeDirectoryMode=2755

UMask=007
PrivateTmp=yes
LimitNOFILE=65535
PrivateDevices=yes
ProtectHome=yes
ReadOnlyDirectories=/
ReadWritePaths=-/var/lib/redis
ReadWritePaths=-/var/log/redis
ReadWritePaths=-/var/run/redis

NoNewPrivileges=true
CapabilityBoundingSet=CAP_SETGID CAP_SETUID CAP_SYS_RESOURCE
MemoryDenyWriteExecute=true
ProtectKernelModules=true
ProtectKernelTunables=true
ProtectControlGroups=true
RestrictRealtime=true
RestrictNamespaces=true
RestrictAddressFamilies=AF_INET AF_INET6 AF_UNIX

# redis-server can write to its own config file when in cluster mode so we
# permit writing there by default. If you are not using this feature, it is
# recommended that you replace the following lines with "ProtectSystem=full".
ProtectSystem=true
ReadWriteDirectories=-/etc/redis

[Install]
WantedBy=multi-user.target
Alias=redis.service

Finally, create a service file in /etc/systemd/system called rqworker@.service to allow ad hoc generation of Redis workers. This file should contain the following lines (replacing <> with appropriate values for your system):

[Unit]
Description=RQ Worker Number %i
After=network.target

[Service]
Type=simple
WorkingDirectory=/home/<username>/<partial path to RoLiM/web directory>
Environment=LANG=en_US.UTF-8
Environment=LC_ALL=en_US.UTF-8
Environment=LC_LANG=en_US.UTF-8
ExecStart=/home/ubuntu/< partial path to virtualenv directory>/patterndetection/bin/<python executable> \
    /home/<username>/<partial path to RoLiM/web/manage.py> \
    rqworker
ExecReload=/bin/kill -s HUP $MAINPID
ExecStop=/bin/kill -s TERM $MAINPID
PrivateTmp=true
Restart=always

[Install]
WantedBy=multi-user.target
```

## Node.js

Install Node.js by following these instructions: https://nodesource.com/blog/installing-node-js-tutorial-ubuntu/

Make sure to replace the version specified in the commands provided in those instructions with the latest version (12.x at the time of writing).

### ReactJS

RoLiM's front end is built in ReactJS. This blog post provides a good overview of using the two together.

https://www.valentinog.com/blog/drf/#django-and-react-together

Follow the instructions provided in that tutorial in order to install React.js and related dependencies.

## Nginx

RoLiM is served using a reverse proxy configuration, where Nginx acts as an HTTP server which listens for and responds to requests from the outside world, and hands them off to Gunicorn via WSGI.

Instructions for installing Nginx on Ubuntu can be found here: https://www.digitalocean.com/community/tutorials/how-to-install-nginx-on-ubuntu-18-04

After installing, navigate to /etc/nginx/sites-available and create a new file called "langelab.org". Add the following lines to that file (replacing <> with appropriate values):

```
server {
  listen 80;
  server_name langelab.org www.langelab.org 206.12.93.194;
  
  client_max_body_size 100M;

  location = /favicon.ico { access_log off; log_not_found off; }
  location /static/ {
    root /home/{username}/<partial path to RoLiM/web/patterndetection/frontend>;
  }

  location / {
    include proxy_params;
    proxy_pass http://unix:/run/gunicorn.sock;
  }
}
```

The IP address on the third line should match the IP address of your server.

Now create a symbolic link to this file in sites-enabled by running the following command:

ln -s langelab.org ../sites-enabled/langelab.org

## Gunicorn

Gunicorn is a WSGI server implemented completely in Python. This resource acts as RoLiM's application server and interacts with Nginx in order to exchange

Create a file at /etc/systemd/system/gunicorn.socket containing the following lines:

```
[Unit]
Description=gunicorn socket

[Socket]
ListenStream=/run/gunicorn.sock

[Install]
WantedBy=sockets.target
```

Now create a file at /etc/systemd/system/gunicorn.service containing the following lines (replacing <> with appropriate values):

```
[Unit]
Description=gunicorn daemon
Requires=gunicorn.socket
After=network.target

[Service]
User=ubuntu
Group=www-data
WorkingDirectory=/home/<username>/<partial path to RoLiM/web>
ExecStart=/home/<username>/<partial path to virtualenv directory>/patterndetection/bin/gunicorn \
    --access-logfile - \
    --workers 3 \
    --bind unix:/run/gunicorn.sock \
    patterndetection.wsgi:application

[Install]
WantedBy=multi-user.target
```

## Systemd

RoLiM's services are managed using systemd. To ensure that the necessary services become available at system start/reboot, run the following commands:

```
sudo systemctl enable rqworker@1.service
sudo systemctl enable gunicorn.socket
sudo systemctl enable gunicorn.service
sudo systemctl enable nginx
```
