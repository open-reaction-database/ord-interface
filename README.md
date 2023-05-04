# ord-interface
Web interface and api for the Open Reaction Database
## Structure
- **/app** - Contains the Vue single page application
  - **/public** - Static SPA assets such as index.html
  - **/src** - Vue source code that gets compiled
    - **/assets** - Images and other static assets
    - **/components** - Reusable Vue components that can be included elsewhere
    - **/router** - Routing files
    - **/styles** - Static .sass files that can be reused throughout the project
    - **/utils** - Static helper .js files
    - **/views** - Routed Vue components, aka "pages"
      - **/browse** - Main browse page of the different reaction sets
      - **/reaction-view** - Display of a single reaction from search results
      - **/search** - Search interface and results
- **/ord_interface** - Contains the Flask app API and legacy interface
  - **/client** - endpoints for the browse and search functionality
  - **/editor** - endpoints for the reaction submission functionality
  - **/visualization** - helper functions for reaction and molecule visuals

## Key Caveats / Constraints
- The Vue front end is dependant on the Flask API. Both must be running for the frontend to work.
- The docker image of the test database must be running on port 5432
- You will need to download the Ketcher interface and extract into appropriate folder (see instructions below) for parts of the interface to work

## How to Deploy
...COMING SOON

## Setup
### 1. Download the repo
```bash
git clone git@github.com:open-reaction-database/ord-interface.git
cd ord-interface
```
### 2. Set up the test database
```bash
cd ./ord_interface
# activate the virtual env of your choice, ex. venv, conda, etc.
# install requirements and run setup script
pip install -e .
./build_test_database.sh 
```
### 3. Set up and run the API via Docker
Note: currently this also runs the old flask ui at port :5001
```bash
# from ./ord_interface
docker build --file Dockerfile -t openreactiondatabase/ord-interface ..
docker compose up
```
  - Leave Docker running in a terminal window.
### 4. Set up and run the Vue SPA
  - Download [Ketcher](https://github.com/epam/ketcher/releases/tag/v2.5.1) (Here's a direct link to the [.zip file](https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip)) and extract the files into `./app/src/ketcher`
  - In a new terminal window:
```bash
cd ./app
# install node packages
npm i 
# run vue spa locally
npm run serve 
```
  - Open [localhost:8080](http://localhost:8080) to view the Vue ORD interface in your browser.
  - The page will reload when you make changes.
  - You may also see any lint errors in the console.