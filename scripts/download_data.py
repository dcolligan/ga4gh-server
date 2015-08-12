"""
Constructs a data source for the ga4gh server by downloading data from
authoritative remote servers.
"""
# TODO
# - refactor into proper inheritence hierarchy
# - need some kind of checkpoint functionality to resume process
#   where it left off since getting a clean run is so rare...
# - implement EBI class; corresponding urls for EBI:
#   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/
#   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import glob
import os

import pysam

import utils


def mkdirAndChdirList(dirNameList):
    """
    For each entry in dirNameList, make a directory if it does not exist and
    cd into it
    """
    for directory in dirNameList:
        mkdirAndChdir(directory)


def mkdirAndChdir(dirName):
    """
    Make a directory if it does not exist and cd into it
    """
    if not os.path.exists(dirName):
        os.mkdir(dirName)
    os.chdir(dirName)


def cleanDir():
    """
    Removes genomic files in current directory
    """
    cwd = os.getcwd()
    utils.log("Cleaning out directory '{}'".format(cwd))
    globs = [
        "*.tbi", "*.vcf", "*.vcf.gz", "*.bam", "*.bam.bai", "*.fasta.gz",
        "*.fasta", "*.fasta.fai"]
    for fileGlob in globs:
        fileNames = glob.glob(fileGlob)
        for fileName in fileNames:
            os.remove(fileName)


class AbstractFileDownloader(object):
    """
    Base class for individual site genome file downloaders
    """
    def __init__(self):
        pass


class NcbiFileDownloader(AbstractFileDownloader):
    """
    Downloads files from NCBI
    """
    def __init__(self, args):
        super(NcbiFileDownloader, self).__init__()
        self.args = args
        self.maxPos = 2**30
        self.minPos = 0
        self.datasetId = 'dataset1'
        self.variantSetId = self.args.source
        self.referenceSetId = 'main'

    def _getFilenames(self):
        baseFileName = (
            "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a"
            ".20130502.genotypes.vcf.gz")
        chrNames = [str(i) for i in range(1, 23)]
        fileNames = [baseFileName.format(chrName) for chrName in chrNames]
        # the X and Y files use a different filename prefix
        fileNames.append(
            'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a'
            '.20130502.genotypes.vcf.gz')
        fileNames.append(
            'ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz')
        return fileNames

    def _getBaseUrl(self):
        baseUrl = (
            "ftp://ftp-trace.ncbi.nih.gov"
            "/1000genomes/ftp/release/20130502/")
        return baseUrl

    def _prepareDir(self):
        dirList = [
            self.args.dir_name, self.datasetId, 'variants',
            self.variantSetId]
        mkdirAndChdirList(dirList)
        cleanDir()

    def _escapeDir(self, levels=4):
        # back to orig dir
        for _ in range(levels):
            os.chdir('..')

    def _downloadIndex(self, fileName):

        @utils.Retry(self.args.num_retries, (IOError))
        def downloadTabixIndex(remoteFileUrl):
            # not sure why an IOError is sometimes thrown here, as if
            # pysam can't find the file
            return pysam.TabixFile(remoteFileUrl)

        @utils.Retry(self.args.num_retries, (IndexError))
        def downloadVariantIndex(remoteFileUrl):
            # not sure why an IndexError is sometimes thrown here
            return pysam.VariantFile(remoteFileUrl)

        baseUrl = self._getBaseUrl()
        remoteFileUrl = os.path.join(baseUrl, fileName)
        utils.log("Downloading index for '{}'".format(remoteFileUrl))
        remoteTabixFile = downloadTabixIndex(remoteFileUrl)
        assert len(remoteTabixFile.contigs) == 1
        # we can't get the header from a TabixFile since that throws an
        # IOError, so we open the same file with the VariantFile interface
        remoteVariantFile = downloadVariantIndex(remoteFileUrl)
        return remoteTabixFile, remoteVariantFile

    def _updatePositions(self, compressedLocalFileName, contig):
        localVariantFile = pysam.VariantFile(compressedLocalFileName)
        localIterator = localVariantFile.fetch(contig)
        for record in localIterator:
            if record.start > self.maxPos:
                self.maxPos = record.start
            if record.start < self.minPos:
                self.minPos = record.start
        localIterator = None
        localVariantFile.close()

    def _processFileName(self, fileName):
        remoteIndexFile = fileName + '.tbi'
        remoteTabixFile, remoteVariantFile = self._downloadIndex(
            fileName)
        contig = remoteTabixFile.contigs[0]
        iterator = remoteTabixFile.fetch(contig)
        localFileName, _ = os.path.splitext(fileName)
        utils.log("Writing '{}'".format(localFileName))
        with open(localFileName, 'w') as localFile:
            headerData = str(remoteVariantFile.header)
            localFile.write(headerData)
            numRecords = 0
            for _ in range(self.args.num_vcf_records):
                data = iterator.next()
                localFile.write(data)
                numRecords += 1
            utils.log("{} records written".format(numRecords))
        iterator = None
        remoteTabixFile.close()
        remoteVariantFile.close()
        os.remove(remoteIndexFile)
        utils.log("Compressing '{}'".format(localFileName))
        utils.runCommand('bgzip -f {}'.format(localFileName))
        compressedLocalFileName = localFileName + '.gz'
        utils.log("Indexing '{}'".format(compressedLocalFileName))
        utils.runCommand('tabix {}'.format(compressedLocalFileName))
        self._updatePositions(compressedLocalFileName, contig)

    def downloadVcfs(self):
        self._prepareDir()
        fileNames = self._getFilenames()
        for fileName in fileNames:
            self._processFileName(fileName)
        self._escapeDir()

    def downloadBams(self):
        dirList = [
            self.args.dir_name, self.datasetId, 'reads', 'low-coverage']
        mkdirAndChdirList(dirList)
        cleanDir()
        studyMap = {
            'HG00096': 'GBR',
            'HG00533': 'CHS',
            'HG00534': 'CHS',
        }
        baseUrl = (
            'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/')
        samples = self.args.samples.split(',')
        for sample in samples:
            samplePath = '{}/alignment/'.format(sample)
            study = studyMap[sample]
            fileName = (
                '{}.mapped.ILLUMINA.bwa.{}.'
                'low_coverage.20120522.bam'.format(sample, study))
            sampleUrl = os.path.join(baseUrl, samplePath, fileName)
            utils.log("Downloading index for '{}'".format(sampleUrl))
            remoteFile = pysam.AlignmentFile(sampleUrl)
            header = remoteFile.header
            utils.log("Writing '{}'".format(fileName))
            localFile = pysam.AlignmentFile(
                fileName, 'wb', header=header)
            self.minPos = 100  # TODO just so this doesn't take forever
            self.maxPos = 200  # TODO remove later
            for reference in remoteFile.references:
                utils.log("reference {}".format(reference))
                iterator = remoteFile.fetch(
                    reference, start=self.minPos, end=self.maxPos)
                numRecords = 0
                for record in iterator:
                    localFile.write(record)
                    numRecords += 1
                utils.log("{} records written".format(numRecords))
            iterator = None
            remoteFile.close()
            localFile.close()
            baiFileName = fileName + '.bai'
            os.remove(baiFileName)
            utils.log("Indexing '{}'".format(fileName))
            # TODO pysam.index(fileName) gives unicode error
            # using command line tool instead
            utils.runCommand("samtools index {}".format(fileName))
        self._escapeDir()

    def downloadRefs(self):
        dirList = [
            self.args.dir_name, 'references', self.referenceSetId]
        mkdirAndChdirList(dirList)
        cleanDir()
        baseUrl = (
            'ftp://ftp-trace.ncbi.nih.gov/'
            '1000genomes/ftp/technical/reference/')
        compressedFastaFileName = 'human_g1k_v37.fasta.gz'
        fastaFileName, _ = os.path.splitext(compressedFastaFileName)
        remoteUrl = os.path.join(baseUrl, compressedFastaFileName)
        maxCoord = 100000000  # TODO change later
        downloader = utils.FastaGunzipFtpFileDownloader(
            remoteUrl, fastaFileName, maxCoord)
        downloader.download()
        utils.log("Indexing {}".format(fastaFileName))
        utils.runCommand("samtools faidx {}".format(fastaFileName))
        self._escapeDir(3)


class EbiFileDownloader(AbstractFileDownloader):
    """
    Downloads files from EBI
    """
    def __init__(self):
        super(EbiFileDownloader, self).__init__()


sources = {
    "ncbi": NcbiFileDownloader,
    "ebi": EbiFileDownloader,
}


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source", default="ncbi", choices=sources.keys(),
        help="the source to download from")
    parser.add_argument(
        "--num-vcf-records", default=5, type=int,
        help="the number of records to pull from each VCF file")
    parser.add_argument(
        "--dir-name", default="ga4gh-downloaded-data",
        help="the name of the directory that the data is downloaded to")
    parser.add_argument(
        "--num-retries", default=3, type=int,
        help="the number of times to retry a failed download")
    parser.add_argument(
        "--samples", default='HG00096,HG00533,HG00534',
        help="a comma-seperated list of samples to download")
    args = parser.parse_args()
    return args


@utils.Timed()
def main(args):
    downloaderClass = sources[args.source]
    downloader = downloaderClass(args)
    # downloader.downloadVcfs()
    # downloader.downloadBams()
    downloader.downloadRefs()


if __name__ == '__main__':
    args = parseArgs()
    main(args)
