export CHROME_BIN=/usr/bin/google-chrome
export DISPLAY=:99.0 
sh -e /etc/init.d/xvfb start
sudo apt-get update
sudo apt-get install -y libappindicator1 fonts-liberation
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
sudo dpkg -i google-chrome*.deb

# chromedriver
wget http://chromedriver.storage.googleapis.com/2.21/chromedriver_linux64.zip
unzip chromedriver_linux64.zip
sudo chmod u+x chromedriver
sudo mv chromedriver /usr/bin/
