gcovr -r . --xml-pretty > report.xml
export CODACY_PROJECT_TOKEN=2280459cadd44135ba6c27aeb3f1d926
bash <(curl -Ls https://coverage.codacy.com/get.sh) report -r report.xml
