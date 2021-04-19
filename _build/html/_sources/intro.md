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

[![badge](https://img.shields.io/badge/MRI%20online%20course-binder-579ACA.svg?logo=data:image/png;base64,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)](https://mybinder.org/v2/gh/zelenkastiot/mri-materials/HEAD)

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