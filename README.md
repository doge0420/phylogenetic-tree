# UPGMA tree builder

# **The script uses [Podman](https://podman.io/) or [Docker](https://www.docker.com/) to run make sure to have one of those engines installed and running before executing the program !**

## How to use ?

First make sure you have inserted either the DNA sequence into the `upgma.py` file (at line 88) or in the distance matrix into the `distances.csv` file. Then see the appropriate category below to run the program. When finished the program will print the intermidiary matrices and output an image with the whole tree in the `/output` folder.

### Run with the script runner (recommended)

#### From DNA sequence

To run the program from raw DNA sequences type and run :

```powershell
python run_upgma.py --mode sequence
```

or simply

```powershell
python run_upgma.py
```

#### From distance matrix

To run the program from a distance matrix type and run :

```powershell
python run_upgma.py --mode matrix
```

Note that you can specify a custom csv file with (default is `distances.csv`) :

```powershell
python run_upgma.py --mode matrix --file [YOUR FILE PATH]
```

**For help run `python run_upgma.py -h`.**

### Run with the raw commands (Not recommended)

#### Build the image with Podman

At the root of the project folder run:

```powershell
podman build -t upgma:latest .
```

#### Run with DNA sequences

```powershell
podman run --rm -v "${PWD}\output:/app/output" upgma:latest
```

#### Run with the distance matrix

Update the `distance.csv` file with the distances and run:

```powershell
podman run --rm -v "${PWD}/output:/app/output" -v "${PWD}/distances.csv:/app/distances.csv" upgma:latest python upgma.py --distances /app/distances.csv
```
