#End 2 End Testing (ngScenario)
To run the end-2-end tests against the application you use Karma or the `runner.html` page.

## Starting the Web Server
In either case you will need the application to be running via the web-server. 
From the root folder of the repository run:

```
node scripts/web-server.js
```
The application should now be available at `http://localhost:3000/app/index.html`


## Testing with Protactor and Selenium server
### Starting the Selenium Server
```
webdriver-manager start --out_dir /Users/jmrodriguez/projects/APPRIS/appris-trunk/server/www/node_modules/protractor/selenium
```
###Start the Protractor test runner using the e2e configuration:
```
protractor config/protractor.conf.js
```

## Testing in a Browser
Browse directly to the e2e test runner page:
```
http://local.es:3000/test/e2e/runner.html
```