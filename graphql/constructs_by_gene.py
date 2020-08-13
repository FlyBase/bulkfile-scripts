#!/usr/bin/env python3
import time
from gql import gql, Client
from gql.transport.requests import RequestsHTTPTransport


def get_client(url, **kwargs):
    """
    Setup the gql transport settings for hitting the FlyBase GraphQL API
    endpoint and return the client.

    :param url: - The URL for the GraphQL endpoint.
    :param kwargs: - Any additional keyword params to pass to the
                     RequestsHTTPTransport class.
    :return: The gql.Client class object.
    """
    flybase_transport = RequestsHTTPTransport(
        url=url,
        use_json=True,
        headers={
            "Content-type": "application/json",
        },
        verify=False,
        retries=3,
        **kwargs
    )
    # Init and return GraphQL client.
    return Client(
        transport=flybase_transport,
        fetch_schema_from_transport=True,
    )


def get_allele_query():
    """
    Function for setting up and returning an allele by gene query.

    This query fetches all construct alleles of a gene specified by the $fbgn parameter.
    Then for each allele it returns the FBal ID, allele symbol, associated constructs,
    regulatory regions, tag uses, and tagged with entities.

    :return: The parsed graphql query object.
    """
    allele_query = '''
    query($fbgn:String!) {
        gene:allelesByGene(fbgn:$fbgn, isConstruct: true) {
            id
            symbol
            alleles {
                id
                symbol
                constructs {
                    id
                    symbol
                }
                regRegions {
                    id
                    symbol
                }
                tagUses {
                    id
                    name
                }
                taggedWith {
                    id
                    symbol
                }
            }
        }
    }
        '''
    # Parse and return the GraphQL query object
    return gql(allele_query)


def main():
    """
    This is a small example script showing how to interface with the FlyBase GraphQL API.
    This script fetches Transgenic construct information for one or more genes. In addition
    the Regulatory region, Tagged with, and Tag Uses are also returned.

    :return: - None
    """

    # List of Genes to process.
    fbgns = ['FBgn0000721', 'FBgn0038032',
             'FBgn0038033', 'FBgn0005626',
             'FBgn0003742', 'FBgn0035697']

    # Get the GQL client and query.
    client = get_client('http://api.flybase.org/graphql')
    query = get_allele_query()

    for fbgn in fbgns:
        # Execute query with fbgn parameter.
        params = {"fbgn": fbgn}
        result = client.execute(query, params)

        # Process all alleles returned for this gene.
        gene = result['gene']
        for allele in gene['alleles']:
            fbal = allele['id']

            # All constructs (FBtp) associated with this allele.
            # This will be a list of dictionaries with 'id' and 'symbol' as keys.
            constructs = allele['constructs']
            print(f"{fbgn} {fbal} Constructs: {constructs}")

            # All regulatory regions (gene; FBgn or tool; FBto) associated with the allele.
            # This will be a list of dictionaries with 'id' and 'symbol' as keys.
            reg_regions = allele['regRegions']

            # All tag uses (controlled vocabulary term; FBcv) associated with the allele.
            # This will be a list of dictionaries with 'id' and 'name' as keys.
            tag_uses = allele['tagUses']

            # All tagged with tools (tool; FBto) associated with the allele.
            # This will be a list of dictionaries with 'id' and 'symbol' as keys.
            tagged_with = allele['taggedWith']

            if len(reg_regions) > 0:
                print(f"{fbgn} {fbal} Reg region: {reg_regions}")
            if len(tag_uses) > 0:
                print(f"{fbgn} {fbal} Tag Uses: {tag_uses}")
            if len(tagged_with) > 0:
                print(f"{fbgn} {fbal} Tagged With: {tagged_with}")

        # Be nice to the GraphQL server.
        time.sleep(1)


if __name__ == '__main__':
    main()
    exit(0)
