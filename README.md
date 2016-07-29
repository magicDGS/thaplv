[![build status](https://api.travis-ci.org/magicDGS/thaplv.svg?branch=master)](https://travis-ci.org/magicDGS/thaplv)
[![Coverage Status](https://coveralls.io/repos/github/magicDGS/thaplv/badge.svg?branch=master)](https://coveralls.io/github/magicDGS/thaplv?branch=master)
[![Language](http://img.shields.io/badge/language-java-brightgreen.svg)](https://www.java.com/)
[![Sputnik](https://sputnik.ci/conf/badge)](https://sputnik.ci/app#/builds/magicDGS/thaplv)
[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


_thaplv_: Tools for Haploid/Haplotype VCFs
=======================================

Tools for Haploid/Haplotype VCFs (_thaplv_) is a software package for analyse haplotype (or haploid) data stored in the
VCF format (see more information about format specifications in genomics [here](http://samtools.github.io/hts-specs/).
The main purpose of this toolkit is include a set of efficient and fast implementations of methods for haplotypes using
the [GATK4 framework](https://github.com/broadinstitute/gatk).

__Note that this project is in early development, so there is no guarantee of maintenance in the current stage.__


## Requirements

* Java 8
* Git
* (Developers) Graddle 2.13 or greater


## Building _thaplv_

For building _thaplv_, run __`./gradlew shadowJar`__. The packaged jar will be in `build/libs/` with the name
`thaplv-<version>-all.jar`.


## Testing _thaplv_

* To run all tests, run __`./gradlew test`__. Test reports will be in `build/reports/tests/index.html`
* To run a subset of tests, use the [test filtering](https://docs.gradle.org/current/userguide/java_plugin.html#test_filtering)
provided by gradle. For example:
  - `./gradlew test --tests *SomeSpecificTest`
  - `./gradlew test --tests all.in.specific.package*`
* To compute coverage reports, run __`./gradlew jacocoTestReport`__. The report will be in `build/reports/jacoco/test/html/index.html`


## Running _thaplv_

For running _thaplv_ after building, run with the command __`java -jar build/libs/thaplv-<version>-all.jar`__.


## Guidelines for developers

* __Pushing directly to master branch is not allowed.__
* It is recomended to name branches with a short form of your name and a explanatory name. Example: dgs_fix_issue30.
* Pull requests should be reviewed before merging by other developer, and rebased/squashed by the author.
* Any new code will require unit/integration tests.
* Use [org.apache.logging.log4j.Logger](https://logging.apache.org/log4j/2.0/log4j-api/apidocs/org/apache/logging/log4j/Logger.html)
for logging.
* Use [TestNG](http://testng.org/doc/index.html) for testing.
* Use [magicDGS Style Guide](https://github.com/magicDGS/styleguide) for code formatting.


## Authors

* Daniel Gomez-Sanchez <daniel.gomez.sanchez@hotmail.es>


## License

Licensed under the [BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause).
See [LICENSE.txt](https://github.com/magicDGS/thaplv/blob/master/LICENSE.txt) file.

The [GATK4 framework](https://github.com/broadinstitute/gatk) is also licensed under the
[BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause).
See their [LICENSE.tx](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT) file for more information.
