#!/usr/bin/env python3
"""Probe W&B for niome Tables / artifacts."""
import json
import urllib.request

GRAPHQL = "https://api.wandb.ai/graphql"


def gql(query: str, variables: dict | None = None) -> dict:
    body = {"query": query}
    if variables:
        body["variables"] = variables
    req = urllib.request.Request(
        GRAPHQL,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=25) as resp:
        return json.loads(resp.read().decode())


def main():
    # Artifacts in project
    q = """
    query {
      project(name: "niome", entityName: "genomes") {
        artifactTypes(first: 10) {
          edges {
            node {
              name
              artifacts(first: 3) {
                edges {
                  node {
                    artifactSequence { name }
                    version
                    createdAt
                    description
                  }
                }
              }
            }
          }
        }
      }
    }
    """
    try:
        out = gql(q)
        print(json.dumps(out, indent=2)[:4000])
    except Exception as e:
        print("artifacts error:", e)

    # Files on latest run
    q2 = """
    query {
      project(name: "niome", entityName: "genomes") {
        run(name: "frf3562w") {
          displayName
          files(names: ["wandb-table.json", "media/table/table.table.json"]) {
            edges { node { name directUrl } }
          }
          outputArtifacts { edges { node { artifactSequenceName version } } }
        }
      }
    }
    """
    try:
        out2 = gql(q2)
        print("\n--- run files ---")
        print(json.dumps(out2, indent=2)[:4000])
    except Exception as e:
        print("files error:", e)

    # Sweep / reports - table URL might be workspace panel
    q3 = """
    query {
      entity(name: "genomes") {
        projects(first: 5) {
          edges {
            node {
              name
              runs(first: 1, order: "-createdAt") {
                edges { node { name displayName } }
              }
            }
          }
        }
      }
    }
    """
    out3 = gql(q3)
    print("\n--- genomes projects ---")
    for e in out3.get("data", {}).get("entity", {}).get("projects", {}).get(
        "edges", []
    ):
        print(" ", e["node"]["name"])


if __name__ == "__main__":
    main()
