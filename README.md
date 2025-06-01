# UPGMA tree builder

## Build the image with Podman

At the root of the project folder run:

```powershell
podman build -t upgma:latest .
```

## Execute the program

### Via DNA sequences

```powershell
podman run --rm -v "${PWD}\output:/app/output" upgma:latest
```

### Via the distance matrix

Update the `distance.csv` file with the distances and run:

```powershell
podman run --rm -v "${PWD}/output:/app/output" -v "${PWD}/distances.csv:/app/distances.csv" upgma:latest python upgma.py --distances /app/distances.csv
```