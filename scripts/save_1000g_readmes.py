"""
Download the AWS 1000genomes READMEs and place them in a folder

Might requre an install of beautiful soup:
    $ pip install beautifulsoup4
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import urlparse

import bs4
import requests


url = 'http://s3.amazonaws.com/1000genomes/'
readmeXmlUrl = url + '?prefix=README'
readmesDirName = '1000genomesREADMEs'


def main():
    if not os.path.exists(readmesDirName):
        os.mkdir(readmesDirName)

    readmeXml = requests.get(readmeXmlUrl)
    soup = bs4.BeautifulSoup(readmeXml.text)
    for contents in soup.find_all('contents'):
        readmeFilename = contents.key.contents[0]
        readmeUrl = urlparse.urljoin(url, readmeFilename)
        readmeRequest = requests.get(readmeUrl)
        readmeFilePath = os.path.join(readmesDirName, readmeFilename)
        print("Writing '{}'".format(readmeFilePath))
        with open(readmeFilePath, 'w') as readmeFile:
            readmeFile.write(readmeRequest.text)


if __name__ == '__main__':
    main()
