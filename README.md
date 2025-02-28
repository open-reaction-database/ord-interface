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
- The Vue front end is dependant on the FastAPI server. Both must be running for the frontend to work.
- The test database must be running on port 5432.
- You will need to download the Ketcher interface and extract into appropriate folder (see instructions below) for 
  parts of the interface to work.

## How to Deploy
...COMING SOON

## Setup

### 1. Download and install

```bash
git clone https://github.com/open-reaction-database/ord-interface
cd ord-interface
# If you are running on Apple silicon, use `conda install postgresql` instead.
conda install -c rdkit rdkit-postgresql
pip install -e '.[tests]'
```

### 2. Set up and run the API

#### Option 1: Docker

```shell
cd ord_interface
./build_test_database.sh
# If you are running on Apple silicon, append `--build-arg="ARCH=aarch_64"` to the next command.
docker build --file Dockerfile -t openreactiondatabase/ord-interface ..
docker compose up
```

#### Option 2: Local

```shell
cd ord_interface/api
ORD_INTERFACE_TESTING=TRUE fastapi dev main.py --port=5000
```

### 3. Set up and run the Vue SPA
  - Download [Ketcher](https://github.com/epam/ketcher/releases/tag/v2.5.1) (Here's a direct link to the [.zip file](https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip)) and extract the files into `./app/src/ketcher`
  - In a new terminal window:

    ```shell
    # you will need node.js and npm installed on your computer first
    cd ./app
    # install node packages
    npm i 
    # run vue spa locally
    npm run serve
    ```

  - Open [localhost:8080](http://localhost:8080) to view the Vue ORD interface in your browser.
  - The page will reload when you make changes.
  - You may also see any lint errors in the console.

## Updating the Schema

### Minor changes such as ENUM additions
  - The Vue app does not require specific modification to show additional ENUM options within existing message fields. We just have to update the package configuration to use the new ord-schema version.
  - Update the ord-schema version number in the [install_requires configuration in ./setup.py](https://github.com/open-reaction-database/ord-interface/blob/aa37f628b176ca241d0701b4df5f6fd7b3079bef/setup.py#L48)
  - Update the [code section in ./ord_interface/Dockerfile](https://github.com/open-reaction-database/ord-interface/blob/aa37f628b176ca241d0701b4df5f6fd7b3079bef/ord_interface/Dockerfile#L57C1-L60C38) which handles the ord-schema install
  - Contact the ORD site administrator to request upload of the new ord-interface version to the staging environment for testing.

### Major changes
  Where there are more extensive changes to the schema, for example addition of new message types, then the Vue app will need to be modified to show these fields. Please discuss any proposed changes with the ORD team.    
