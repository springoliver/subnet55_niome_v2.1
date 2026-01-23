#!/usr/bin/env python3
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
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def main():
    run = "frf3562w"
    q = """
    query($entity: String!, $project: String!, $run: String!) {
      project(name: $project, entityName: $entity) {
        run(name: $run) {
          displayName
          historyKeys
          sampledHistory(specs: [{maxSamples: 20}]) {
            keys
            values
          }
        }
      }
    }
    """
    out = gql(q, {"entity": "genomes", "project": "niome", "run": run})
    run_obj = out["data"]["project"]["run"]
    print("displayName:", run_obj["displayName"])
    print("historyKeys:", json.dumps(run_obj.get("historyKeys"), indent=2)[:2000])
    sh = run_obj.get("sampledHistory") or []
    for i, block in enumerate(sh):
        print(f"\n--- sampledHistory block {i} ---")
        print("keys:", block.get("keys"))
        for v in (block.get("values") or [])[:3]:
            s = json.dumps(v) if not isinstance(v, str) else v
            print(s[:800])


if __name__ == "__main__":
    main()
