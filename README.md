# UPGMA tree builder

# Disclaimer

This program was developed and tested on **Windows 11 Pro** with **Python 3.13.4**. If you're using a different operating system or Python version, you may encounter compatibility issues.

To minimize such problems, the `upgma.py` script is designed to run inside a **containerized environment** using either Docker or Podman. This ensures a consistent runtime setup across platforms.

Please note that while the script has been tested on several example cases and produced correct results, I **cannot guarantee** its correctness in all scenarios. Use this tool **at your own risk** and validate results independently if accuracy is critical.

## ⚠️ Prerequisites 

- A somewhat modern version of python installed (>= 3.10 to be safe).
- Either [Podman](https://podman.io/) or [Docker](https://www.docker.com/) installed and running in the background.

## How to use ?

First make sure you have inserted either the DNA sequence into the `upgma.py` file (at line 88) or in the distance matrix into the `distances.csv` file. Then see the appropriate category below to run the program. When finished, the program will print the intermediary matrices and output an image with the whole tree in the `/output` folder.

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
