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

## How to Deploy
...COMING SOON

## Setup
1. Download the repo
```shell
git clone git@github.com:open-reaction-database/ord-interface.git
cd ord-interface
```
2. Set up the test database
```shell
cd ./ord_interface
# activate the virtual env of your choice
pip install -e .
./build_test_database.sh
docker run -d -p 5432:5432 openreactiondatabase/ord-postgres:test #creates docker image of test db
```
3. Set up and run the Flask API (Don't forget to use the virtual environment)
```shell
# from ./ord_interface
POSTGRES_USER=postgres POSTGRES_PASSWORD=postgres FLASK_APP=interface.py FLASK_ENV=development python -m flask run
```
  - Leave the Flask app running in a terminal window.
4. Set up and run the Vue SPA
  - Download [https://github.com/epam/ketcher/releases/tag/v2.5.1](Ketcher) (Here's a direct link to the [https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip](.zip file)) and extract the files into `./app/src/ketcher`
  - In a new terminal window:
```shell
cd ./app
npm i #install node packages
npm run serve #runs vue spa locally
```
  - Open [http://localhost:8080](http://localhost:8080) to view the Vue ORD interface in your browser.
  - The page will reload when you make changes.
  - You may also see any lint errors in the console.