#!/usr/bin/env python3
"""Fetch W&B niome validator run metrics (public GraphQL)."""
import json
import urllib.request

GRAPHQL = "https://api.wandb.ai/graphql"
ENTITY = "genomes"
PROJECT = "niome"
# Validator from Results/miner.json
RUN_NAME = "f14tpkye"  # validator-119-2.1.0


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
    # Run summary + history keys
    q1 = """
    query($entity: String!, $project: String!, $run: String!) {
      project(name: $project, entityName: $entity) {
        run(name: $run) {
          name
          displayName
          createdAt
          summaryMetrics
          config
        }
      }
    }
    """
    out = gql(
        q1,
        {"entity": ENTITY, "project": PROJECT, "run": RUN_NAME},
    )
    run = out.get("data", {}).get("project", {}).get("run") or {}
    print("=== Run", run.get("displayName"), RUN_NAME, "===")
    sm = run.get("summaryMetrics")
    if sm:
        if isinstance(sm, str):
            sm = json.loads(sm)
        print("summaryMetrics keys:", sorted(sm.keys())[:40])
        for k in sorted(sm.keys()):
            if any(
                x in k.lower()
                for x in (
                    "score",
                    "uid",
                    "miner",
                    "vcf",
                    "final",
                    "task",
                    "variant",
                    "top",
                )
            ):
                print(f"  {k}: {sm[k]}")

    # History samples
    q2 = """
    query($entity: String!, $project: String!, $run: String!) {
      project(name: $project, entityName: $entity) {
        run(name: $run) {
          historyKeys
        }
      }
    }
    """
    try:
        out2 = gql(q2, {"entity": ENTITY, "project": PROJECT, "run": RUN_NAME})
        keys = out2.get("data", {}).get("project", {}).get("run", {}).get(
            "historyKeys"
        )
        if keys:
            print("\nhistoryKeys (score-related):")
            for k in sorted(keys):
                if any(
                    x in k.lower()
                    for x in ("score", "uid", "miner", "vcf", "final", "task")
                ):
                    print(f"  {k}")
    except Exception as e:
        print("historyKeys error:", e)

    # Latest runs list
    q3 = """
    {
      project(name: "niome", entityName: "genomes") {
        runs(first: 8, order: "-createdAt") {
          edges {
            node {
              name
              displayName
              createdAt
            }
          }
        }
      }
    }
    """
    out3 = gql(q3)
    print("\n=== Recent validator runs ===")
    for edge in out3["data"]["project"]["runs"]["edges"]:
        n = edge["node"]
        print(f"  {n['createdAt']}  {n['displayName']}  ({n['name']})")


if __name__ == "__main__":
    main()
