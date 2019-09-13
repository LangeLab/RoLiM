import React from "react";
import ReactDOM from "react-dom";
import DataProvider from "./DataProvider";
import Table from "./Table";
import Form from "./Form";
import JobForm from "./JobForm";
import JobStatus from "./JobStatus";

const App = () => (
    <React.Fragment>
      {/*<DataProvider endpoint="api/job/" 
                render={data => <Table data={data} />} />*/}
      <JobStatus endpoint="api/job/" />
      <JobForm endpoint="api/job/" />
    </React.Fragment>
);
const wrapper = document.getElementById("app");
wrapper ? ReactDOM.render(<App />, wrapper) : null;