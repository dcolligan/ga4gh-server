"""
An end to end test which tests:
- client cmd line parsing
- client operation
- client logging
- server cmd line parsing
- server operation
- simulated variantSet backend
- server logging
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import server_test
import client

import ga4gh.datamodel as datamodel
import tests.paths as paths


class TestGestalt(server_test.ServerTest):
    """
    An end-to-end test of the client and server
    """
    def testEndToEnd(self):
        self.simulatedDatasetId = datamodel.DatasetCompoundId.obfuscate(
            paths.simulatedDatasetId)
        self.simulatedVariantSetId = \
            datamodel.VariantSetCompoundId.obfuscate(
                paths.simulatedVariantSetId)
        self.simulatedReadGroupId = datamodel.ReadGroupCompoundId.obfuscate(
            paths.simulatedReadGroupId)
        self.simulatedReferenceSetId = \
            datamodel.ReferenceSetCompoundId.obfuscate(
                paths.simulatedReferenceSetId)
        self.simulatedReferenceId = \
            datamodel.ReferenceCompoundId.obfuscate(
                paths.simulatedReferenceId)
        self.simulatedVariantAnnotationSetId = \
            datamodel.VariantAnnotationSetCompoundId.obfuscate(
                paths.simulatedVariantAnnotationSetId)
        self.client = client.ClientForTesting(self.server.getUrl())
        self.runVariantsRequest()
        self.assertLogsWritten()
        self.runReadsRequest()
        self.runReferencesRequest()
        self.runVariantSetsRequestDatasetTwo()
        self.client.cleanup()

    def assertLogsWritten(self):
        serverOutLines = self.server.getOutLines()
        serverErrLines = self.server.getErrLines()
        clientOutLines = self.client.getOutLines()
        clientErrLines = self.client.getErrLines()

        # nothing should be written to server stdout
        self.assertEqual(
            [], serverOutLines,
            "Server stdout log not empty")

        # server stderr should log at least one response success
        responseFound = False
        for line in serverErrLines:
            if ' 200 ' in line:
                responseFound = True
                break
        self.assertTrue(
            responseFound,
            "No successful server response logged to stderr")

        # client stdout should not be empty
        self.assertNotEqual(
            [], clientOutLines,
            "Client stdout log is empty")

        # number of variants to expect
        expectedNumClientOutLines = 2
        self.assertEqual(len(clientOutLines), expectedNumClientOutLines)

        # client stderr should log at least one post
        requestFound = False
        for line in clientErrLines:
            if 'POST' in line:
                requestFound = True
                break
        self.assertTrue(
            requestFound,
            "No request logged from the client to stderr")

    def runVariantsRequest(self):
        self.runClientCmd(
            self.client,
            "variants-search",
            "-s 0 -e 2 -V {}".format(self.simulatedVariantSetId))

    def runVariantAnnotationsRequest(self):
        self.runClientCmd(
            self.client,
            "variantannotations-search",
            "{}".format(self.simulatedVariantAnnotationSetId))

    def runReadsRequest(self):
        args = "--readGroupIds {} --referenceId {}".format(
            self.simulatedReadGroupId, self.simulatedReferenceId)
        self.runClientCmd(self.client, "reads-search", args)

    def runReferencesRequest(self):
        referenceSetId = self.simulatedReferenceSetId
        referenceId = self.simulatedReferenceId
        cmd = "referencesets-search"
        self.runClientCmd(self.client, cmd)
        cmd = "references-search"
        args = "--referenceSetId={}".format(referenceSetId)
        self.runClientCmd(self.client, cmd, args)
        cmd = "referencesets-get"
        args = "{}".format(referenceSetId)
        self.runClientCmd(self.client, cmd, args)
        cmd = "references-get"
        args = "{}".format(referenceId)
        self.runClientCmd(self.client, cmd, args)
        cmd = "references-list-bases"
        args = "{}".format(referenceId)
        self.runClientCmd(self.client, cmd, args)

    def runVariantSetsRequestDatasetTwo(self):
        cmd = "variantsets-search"
        args = "--datasetId {}".format(self.simulatedDatasetId)
        self.runClientCmd(self.client, cmd, args)
