$(document).ready(function(){
    $('#resultTable').DataTable({
        "bPaginate": false,
         columnDefs: [
            { type: 'natural', targets: '_all' }
         ]
    });
});

$(document).ready(function(){
  $( "#invert_button" ).click(function() {
    console.log("INVERT!")
    // the css we are going to inject
    var css = 'html {-webkit-filter: invert(100%);background-color: black;' +
        '-moz-filter: invert(100%);background-color: black;' +
        '-o-filter: invert(100%);background-color: black;' +
        '-ms-filter: invert(100%);background-color: black; }',

    head = document.getElementsByTagName('body')[0],
    style = document.createElement('style');

    // a hack, so you can "invert back" clicking the bookmarklet again
    if (!window.counter) { window.counter = 1;} else  { window.counter ++;
    if (window.counter % 2 == 0) { var css ='html {-webkit-filter: invert(0%);background-color: white; -moz-filter:    invert(0%);background-color: white; -o-filter: invert(0%);background-color: white; -ms-filter: invert(0%);background-color: white; }'}
     };

    style.type = 'text/css';
    if (style.styleSheet){
    style.styleSheet.cssText = css;
    } else {
    style.appendChild(document.createTextNode(css));
    }

    //injecting the css to the head
    head.appendChild(style);
  });
});