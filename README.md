# ASOstool-v2 — User & Developer Guide

## Setup guide

### 👥 Regular user setup (Shiny server only)

This section is for end-users who only need to run the Shiny application. No RStudio or Git is required. 

#### 1. Build the Docker image 

```docker build -t asostoolv2 .```

#### 2. Start the container

```         
docker run -d \
    -p 3838:3838 \
    -v $(pwd):/srv/shiny-server/ASOstool-v2 \
    asostoolv2
```

#### 3. Open the Shiny App

Open:

👉 <http://localhost:3838/ASOstool-v2/Simpel/> 

You're done!

#### 4. Stop the container

Stop the container after you're done to free up ports.
List running containers:

```
docker ps
```

Stop it:

```
docker stop <container_id>
```

Remove it (optional):

```
docker rm <container_id>
```

### 💻 Developer setup (RStudio + Git + SSH)

This section is for developers contributing to the project.

#### 1. Create a persistent volume for RStudio home

```
docker volume create rstudio-home
```

#### 2. Start the development container

Mount:

-   Project folder → RStudio
-   Project folder → Shiny
-   Persistent RStudio home folder

```         
docker run -d \
  -p 8787:8787 \
  -p 3838:3838 \
  -v $(pwd):/home/rstudio/ASOstool-v2 \
  -v $(pwd):/srv/shiny-server/ASOstool-v2 \
  -v rstudio-home:/home/rstudio \
  asostoolv2
```

#### 3. Log into RStudio

Open:

👉 <http://localhost:8787>

Credentials:

-   Username: rstudio
-   Password: rstudio

#### 4. Set up SSH keys (first time only)

Open the Terminal inside RStudio.

##### 4.1 Start SSH agent

```eval "$(ssh-agent -s)"```

##### 4.2 Generate SSH key

Replace example email with github email

```ssh-keygen -t ed25519 -C "your_email@example.com"```

Press ENTER for all prompts.

##### 4.3 Add the key to ssh-agent

```ssh-add ~/.ssh/id_ed25519```

##### 4.4 Add your public key to GitHub

```cat ~/.ssh/id_ed25519.pub```

Copy → GitHub →

**Settings** → **SSH and GPG Keys** → **New SSH key**

#### 5. Configure Git user (first time only)

```         
git config --global user.name "Your Name"`
git config --global user.email "you@example.com"
```

These persist thanks to the rstudio-home volume.

#### 6. Fix Push/Pull — convert HTTPS to SSH

If your remote url is set trough https you need to switch to SSH

```
cd /home/rstudio/ASOstool-v2
git remote set-url origin git@github.com:JimReijnhoudt/ASOstool-v2.git
```

Check:

```
git remote -v
```

You should now see:

```
origin  git@github.com:JimReijnhoudt/ASOstool-v2.git (fetch)
origin  git@github.com:JimReijnhoudt/ASOstool-v2.git (push)
```

Now the Push and Pull buttons in RStudio Git tab become active.

#### 7. Fix: “All files show as modified”

Run once:

```
git config core.autocrlf input
```

#### 8. Stop the development container

Stop the container after you're done to free up ports.
List running containers:

```
docker ps
```

Stop it:

```
docker stop <container_id>
```

Remove it (optional):

```
docker rm <container_id>
```

### ✔ Setup Summary Table

| Action                                                     | Required Once | Required Each Restart |
| ----------------------------------------------------------- | ------------- | --------------------- |
| Create `rstudio-home` volume                                | ✔             | ❌                     |
| Generate SSH key                                            | ✔             | ❌                     |
| Add SSH key to GitHub                                       | ✔             | ❌                     |
| Set Git username/email                                      | ✔             | ❌                     |
| Convert remote to SSH                                       | ✔             | ❌                     |
| Start ssh-agent + ssh-add (If using password for SSH key)   | ❌             | ✔                     |
| Run docker container                                        | ❌             | ✔                     |
| Stop container                                              | ❌             | ✔                     |
