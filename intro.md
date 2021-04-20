# Introduction


Welcome to an interactive **MRI course**! This online course is divided in 4 chapters (for more info see 
<font> <b> 📕 Course chapters </b> </font>) that demonstrate different MRI techniques and display the results
in *Jupyer notebooks*. The notebooks are written in *Julia* (1.4.1) and *Python* (3.7).  

If you are curious on how to use this *Jupyter-book* interface, download the codes or reproduce the figures see these sections:

<details><summary><font> <b> 📕 Course chapters </b> </font> </summary><br>

**1: Sensitivity encoded MRI reconstruction**
- Preview: https://nbviewer.jupyter.org/github/nankueichen/education_SENSE_MRI/blob/main/main.ipynb
- Complete files: https://github.com/nankueichen/education_SENSE_MRI

**2: Magnitude and phase data processing for multi-TE gradient-echo MRI**
- Preview: https://nbviewer.jupyter.org/github/nankueichen/multi_TE_MRI/blob/main/main.ipynb
- Complete files: https://github.com/nankueichen/multi_TE_MRI

**3: DTI data processing**
- Preview: https://nbviewer.jupyter.org/github/nankueichen/education_DTI_processing/blob/main/main.ipynb
- Complete files: https://github.com/nankueichen/education_DTI_processing

**4: RF pulse design:**
- Preview: https://nbviewer.jupyter.org/github/sequintoa/BME639-RF/blob/master/Lab/RFPulse.ipynb
- Complete files: https://github.com/sequintoa/BME639-RF/tree/master/Lab

---
</details>

<details><summary><font><b>🐳 Docker enviroment</b> </font></summary><br>

Dockerfile for running a Docker image able to run [SoS](https://vatlab.github.io/sos-docs/) Jupyter notebooks (**Julia**: 1.4.1, **Python**: 3.7) 

**Run Docker locally** 


If you have Docker installed on your computer and running, you can run the code 
in the same environment described in this repository using `repo2docker`. 

1. Simply install `repo2docker` from pyPI: 
```
pip install jupyter-repo2docker
```
2. Run the following command in your terminal:
```
jupyter-repo2docker https://github.com/neurolibre/myelin-meta-analysis
```

After building (it might take a while!), it should output in your terminal 
something like:

```
Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://0.0.0.0:36511/?token=f94f8fabb92e22f5bfab116c382b4707fc2cade56ad1ace0
```

This should start a Jupyter session on your browser and make all the resources 
you see when you [launch a Binder](https://mybinder.org/v2/gh/neurolibre/myelin-meta-analysis/master) for this repository. 

To re-use your container built by repo2docker, do the following: 

1. Run `docker images` command and copy the `IMAGE ID` to your clipboard 
2. Run the following command to start the container:
```
docker run -it --rm -p 8888:8888 `PASTE IMAGE ID HERE` jupyter notebook --ip 0.0.0.0
```
---
</details>

<details><summary><font><b>☁️ Run on the cloud </b> </font> </summary><br>

You can use <code> Live Code </code> or <code>Launch in Binder</code> buttons 
at the top of each page of the <a href="">Jupyter Book</a>.

Alternatively, you can start a Binder session by clicking the badge below: 

[![badge](https://raw.githubusercontent.com/Notebook-Factory/mri-course/8e3a28f4f0c9c32ee2f2b9df0055271ebcb91a9b/images/MRI-online.svg)](https://mybinder.org/v2/gh/zelenkastiot/mri-materials/HEAD)

---
</details>

<details><summary><font><b>▶️YouTube guide </b> </font> </summary><br>

Insert video...

--- 
</details>

<br>

<p align="center">
<img src="https://avatars3.githubusercontent.com/u/63861117?s=200&v=4" style="width:34px;"></img> <br>
This repository is created by <a href="https://github.com/Notebook-Factory">Notebok-Factory</a>. 
</p>