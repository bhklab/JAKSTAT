const path = require('path');
const fs = require('fs');

(async() => {
    let dirs = fs.readdirSync("/Users/minorunakano/Documents/JAKSTAT/tcl38/output");
    dirs = dirs.filter(dir => dir !== ".DS_Store")
    for(const dir of dirs){
        let runInfo = fs.readFileSync(`/Users/minorunakano/Documents/JAKSTAT/tcl38/output/${dir}/run_info.json`);
        runInfo = JSON.parse(runInfo);
        console.log(`${dir} ${runInfo.p_pseudoaligned} ${runInfo.p_pseudoaligned < 20.0 ? 'Low!' : ''}`);
    }
})()