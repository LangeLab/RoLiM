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

Python should be installed as an altinstall rather than replacing the default system Python version in order to prevent possible system-wide compatibility issues.

## MariaDB

RoLiM uses the MariaDB implementation of MySQL. Installation instructions can be found here: https://www.digitalocean.com/community/tutorials/how-to-install-mariadb-on-ubuntu-20-04

### Configuration

After installing, start a MySQL session by running "sudo mysql".

Next, create a database called "patterndetection" by running the following command:


Create a new user named "patterndetection" by running the following command:



Create a second database called "merops" by running the following command:

Now create a new user "zigzag" with full privileges on the merops database by running the following commands:

CREATE USER "zigzag"@"localhost" identified by "zigzagpass";



### Installation notes

Python's MySQL connector requires the C client library API which can be installed by running "sudo apt-get install libmysqlclient-dev"

## Redis

https://www.digitalocean.com/community/tutorials/how-to-install-and-secure-redis-on-ubuntu-20-04

Require password

## Node.js

Install Node.js by following these instructions: https://nodesource.com/blog/installing-node-js-tutorial-ubuntu/

Replace the version specified in the commands provided in those instructions with the latest version (12.x at the time of writing).

### ReactJS

RoLiM's front end is built in ReactJS. This blog post provides a good overview of using the two together.

https://www.valentinog.com/blog/drf/#django-and-react-together


Follow the instructions provided in that tutorial in order to install React.js and related dependencies.

## Nginx

RoLiM is served using a reverse proxy configuration, where Nginx acts as an HTTP server which listens for and responds to requests from the outside world, and hands them off to Gunicorn via WSGI.

Instructions for installing Nginx on Ubuntu can be found here: https://www.digitalocean.com/community/tutorials/how-to-install-nginx-on-ubuntu-18-04

## Gunicorn

Gunicorn is a WSGI server implemented completely in Python. This resource acts as RoLiM's application server, and interacts with Nginx in order to exchan