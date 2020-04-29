import React from "react";
import ReactDOM from "react-dom";
import DataProvider from "./DataProvider";
import Table from "./Table";
import JobForm from "./JobForm";

const App = () => (
  <div className="container">
  	<div>
  		<h1 class="title">RoLiM</h1>
  		<h2 class="st">Robust detection of linear motifs in sequence data.</h2>
  	</div>
  	<br />
  	<JobForm endpoint="/api/job/" />
  </div>
);
const wrapper = document.getElementById("app");
wrapper ? ReactDOM.render(<App />, wrapper) : null;
