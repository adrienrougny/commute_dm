import json

with open("results_annotation_partly_reviewed.json") as f:
    data = json.load(f)
    print(len(data))
    data_filtered = [_ for _ in data if _["comments"]]
    print(len(data_filtered))

with open("results_annotation_filtered.json", "w") as f:
    json.dump(data_filtered, f, indent=2)

for _ in data:
    if not _["comments"] and not _["additional_annotations"]:
        print(_)
