# *BAGS: Barcode, Audit & Grade System*

<div style="border: 2px solid red; padding: 15px; background-color: #ffe6e6; font-size: 18px;">
  <h1 style="color:red; text-align:center;">⚠️ Important Notice ⚠️</h1>

  <p>
    <strong>BAGS is currently not operational.</strong>
  </p>

  <p>
    BAGS was originally developed to retrieve barcode records from BOLD Systems using the BOLD v4 API and the R package <code>bold</code>. Following the transition of BOLD Systems to version 5 and its new API infrastructure, the v4 services on which BAGS relies no longer provide the data access required by the application. Consequently, BAGS is currently unable to retrieve the BOLD records needed to perform its analyses and grade assignments.
  </p>

  <p>
    In addition, the R package <code>bold</code>, which BAGS uses for data retrieval, was removed from the CRAN repository in August 2024. Although archived versions of the package remain available, they were developed for the BOLD v4 API and therefore do not restore the current functionality of BAGS.
  </p>

  <p>
    Users interested in applying a similar barcode reference-library auditing approach using the current BOLD infrastructure may also explore <strong>BOLDcuratoR</strong>, a tool developed within the <strong>Biodiversity Genomics Europe (BGE)</strong> project that includes the BAGS A–E grade assignment. The tool and its source code are available at
    <a href="https://github.com/bge-barcoding/BOLDcuratoR" target="_blank">https://github.com/bge-barcoding/BOLDcuratoR</a>.
  </p>

  <p>
    <strong>We are currently exploring the feasibility of adapting BAGS to the BOLD v5 API.</strong> This would require updating the data-retrieval workflow and other parts of the application to accommodate the structure and endpoints of the new API. This page will be updated if a compatible version becomes available.
  </p>
</div>
