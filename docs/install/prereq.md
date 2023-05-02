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
    ![x](../../docs/images/assets/git_ssh_key.png)
    - Click `Add SSH Key`.
  - You're all set for now with adding ssh keys to github.


## AWS Console Stuff

### Regions
  - Please create all resources in the us-west-2 (Oregon) region.  This may be changed, but all instructions will assume this region. If you choose a different region, make sure all EC2 instance types used in the daylily ephemeral cluster (DEC) config are also available... or edit out the ones which are not (but it will be easier to stick with us-west-2 for now).

### Key Pair
- Daylily will locally utilize the [aws cli](https://docs.aws.amazon.com/cli/index.html) and [aws parallelcluster cli](https://docs.aws.amazon.com/parallelcluster/latest/ug/commands-v3.html) to create a DEC. 
  > **Note**
  > These credentials are not used beyond the user account you use to configure and launch DEC's, IAM roles are used for all other aspects of AWS interactions, these credentials are not copied elsewhere.

-  [Instructions for creating a key pair](https://docs.aws.amazon.com/parallelcluster/latest/ug/set-up-keypair.html):
    1. Go to the AWS Management Console and log in to your account.
    2. Click on your username on the top right corner and select "My Security Credentials" from the dropdown menu.
    3. On the left-hand side of the page, click on "Access keys (access key ID and secret access key)".
    4. Click on the "Create New Access Key" button.
    5. In the pop-up window that appears, click on the "Download .csv file" button to download a file containing your new access key ID and secret access key. Make sure to keep this file in a secure location as it cannot be retrieved later.
    6. After downloading the .csv file, click on "Close" to return to the Access keys page.
    7. Your new access key and secret access key will be displayed on the page. You can also copy them from this page if needed. 
    > **Note**
    > You only have one chance to save these keys. If you do not do so, or loose them, you will need to create a new keys.
- Save your `access key` and `secret access key` for latter, the aws cli will ask for them when you first are prompted to use it.

### PEM File
  - As with the AWS credential keys, a `PEM` file will be needed to create and access the headnode of the DEC. 
    > **Note**
    > Theis `PEM` file is not used beyond the user account you use to configure and launch DEC's, IAM roles are used for all other aspects of AWS interactions, this `PEM` file is not copied elsewhere.
  - [Create and download](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#having-ec2-create-your-key-pair) to your `PEM` to your `~/.ssh/` directory:
    1. Log in to the AWS Management Console and go to the EC2 dashboard.
    2. Click on the "Key Pairs" link in the left-hand navigation pane.
    3. Click on the "Create Key Pair" button.
    4. In the "Create Key Pair" dialog box, enter a name for your key pair in the "Key pair name" field. Make sure to choose a descriptive name that you will remember later.
    5. Choose the "PEM" format from the "File format" dropdown menu.
    6. Click on the "Create Key Pair" button.
    7. Your key pair will be created, and a private key file with a .pem file extension will be downloaded to your local machine. Save this file in a secure location as you will need it to access your instances.
    8. Download and save the `PEM` file to your `~/.ssh/` directory. *It should be named PEMKEYNAME.pem*.
    > **Note**
    > You only have one chance to download this key. If you do not do so, or loose the file, you will need to create a new `PEM`.

