import sys
import requests
import yaml
from IPython import embed

# Load the Conda YAML file
with open(sys.argv[1], "r") as f:
    conda_yaml = yaml.safe_load(f)

# Extract the names of the open source tools
packages = conda_yaml.get("dependencies", [])
package_names = [p.split("=")[0] for p in packages if isinstance(p, str)]

# Query the Crossref API for scientific citations
citations = []
for package_name in package_names:
    query = f'title:"{package_name}"'
    response = requests.get(f"https://api.crossref.org/works?query.container-title={query}&rows=1")
    if response.status_code == 200:
        results = response.json().get("message", {}).get("items", [])
        for result in results:

            try:
                print(result)
                doi = result.get("DOI", "")
                
                try:
                    title = result.get("container-title", [])[0]
                except Exception as e1:
                    print('E1', e1)
                    try:
                        title = result.get("title", [])[0]
                    except Exception as e2:
                         print('E2', e2)
                         title='NA'
                author = ", ".join(a.get("given", "") + " " + a.get("family", "") for a in result.get("author", []))
                citation = f"{author} ({doi}): {title}"
                citations.append(citation)
            except Exception as e:
                print('E0', e)
                embed()
                raise

# Print the citations
for citation in citations:
    print(citation)
