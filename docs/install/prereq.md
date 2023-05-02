# Daylily Pre-requisites

## Github Account & Stored Local id_rsa.pub
  - Create you local ssh key if you have not done so already.
  ```bash
ssh-keygen 
# Generating public/private rsa key pair.
Enter file in which to save the key (/Users/day/.ssh/id_rsa): 
# Created directory '/Users/day/.ssh'.
Enter passphrase (empty for no passphrase): 
Enter same passphrase again: 
# Your identification has been saved in /Users/day/.ssh/id_rsa
# Your public key has been saved in /Users/day/.ssh/id_rsa.pub
  ```
  
  - Capture Your id_rsa.pub Key.
  ```bash
  head ~/.ssh/id_rsa.pub
  # << your key displayed here >>
  # Double click the key and copy to clipboard. If you use head, it should copy w/out linebreaks.
  ```
  
  - Save Key To github.
    - After logging in to [github](www.github.com), go to the [settings](https://github.com/settings/profile)->[ssh and gpg keys](https://github.com/settings/keys).
    - Click `New SSH Key`.
    - Give the key a unique name, and paste your id_rsa.pub key into the `Key` textarea.
    ![](docs/images/assets/git_ssh_key.png)
    - Click `Add SSH Key`.
  - You're all set for now with adding ssh keys to github.


## 
