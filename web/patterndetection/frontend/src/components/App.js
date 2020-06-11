import React from "react";
import ReactDOM from "react-dom";
import DataProvider from "./DataProvider";
import Table from "./Table";
import JobForm from "./JobForm";

const App = () => (
  <div className="container">
  	<header>
  		<h1 class="title"><a href="https://langelab.med.ubc.ca/">Lange Lab</a></h1>
  		<h2 class="title">RoLiM</h2>
  		<h3 class="subtitle">Robust detection of linear motifs in sequence data.</h3>
  	</header>
  	<br />
  	<JobForm endpoint="/api/job/" />
  </div>
);
const wrapper = document.getElementById("app");
wrapper ? ReactDOM.render(<App />, wrapper) : null;
