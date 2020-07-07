/*
 * Import libraries
 */


/*
 * Global variables
 */
// let wk_file = "../data/FunctionalTableLong.xlsx"

// convert text variable separating by delimiter to array of array
function convert_txt_to_json(txt, sep="\t") {
    let report = [];
    let lines = txt.split("\n");
    for (let i=0; i<lines.length; i++) {
        let cols = lines[i].split(sep);
        let rep = [];
        for (let j=0; j<cols.length; j++) {
            let col = cols[j];
            rep.push(col);
        }
        if ( rep.length >=1) report.push(rep);
    }
    return report;
}

// create handsontable from given data
function create_table(table_id, data) {
  var example1 = document.getElementById(table_id);
  
  var hot = new Handsontable(example1, {
    data: data.slice(1),
    colHeaders: data[0],
    rowHeaders: true,
    filters: true,
    dropdownMenu: ['filter_by_value', 'filter_action_bar'],
    columnSorting: true,
    manualColumnResize: true,
    minSpareRows: 0,
    readOnly: true, // make table cells read-only
    contextMenu: false, // disable context menu to change things
    // disableVisualSelection: true, // prevent user from visually selecting
    // manualColumnResize: false, // prevent dragging to resize columns
    // manualRowResize: false, // prevent dragging to resize rows
    comments: false, // prevent editing of comments
    licenseKey: 'non-commercial-and-evaluation'    
  });
}